#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

struct Pixel {
    int r;
    int g;
    int b;

    Pixel& operator+=(const Pixel& other){
        r += other.r;
        g += other.g;
        b += other.b;

        return *this;
    }

    Pixel& operator/=(int value){
        r /= value;
        g /= value;
        b /= value;

        return *this;
    }

};

struct PixelSum {
    int count;
    Pixel pixel;
};

using Pixels = std::vector<Pixel>;
using PixelsSums = std::vector<PixelSum>;

struct KMeansParams {
    int start;
    int end;
    int centroids;
    int iterations;
    Pixels image;
};

struct KMeansParamsFormMainProcess {
    int centroids;
    int iterations;
    int processes;
    Pixels image;
    Pixels centroidsPoints;
};

struct KMeansMpiParams {
    MPI_Datatype pixelDt;
    MPI_Datatype sumPixelDt;
    int tag;
};

struct MpiReceiveParameters {
    MPI_Datatype type;
    int tag;
    MPI_Status* status;
};

Pixels getImageFromFile(std::string fileName, int width, int hight){
    std::fstream file{fileName};
    Pixels  output;
    for (int i = 0; i < width * hight; ++i){
        int r, g, b;
        file >> r >> g >> b;
        output.push_back({r, g, b});
    }
    return output;
}

void imageToFile(std::string fileName, const Pixels& image){
    std::fstream file{fileName};
    for (int i = 0; i < image.size(); ++i){
        file << image[i].r << ' ' << image[i].g << ' ' << image[i].b << ' ';
    }
}

double euclideanDistance(Pixel imagePixel, Pixel clusterPixel){
    return std::sqrt(std::pow(imagePixel.r - clusterPixel.r, 2) +
                     std::pow(imagePixel.g - clusterPixel.g, 2) +
                     std::pow(imagePixel.b - clusterPixel.b, 2));
}

int assignCluster(Pixel imagePixel, const Pixels& clusters)
{
    std::vector<int> distances;
	for (const auto& clusterPixel : clusters)
	{
		int distance = euclideanDistance(imagePixel, clusterPixel);
        distances.push_back(distance);

	}
	return std::distance(std::begin(distances), std::min_element(std::begin(distances), std::end(distances)));
}

Pixels getRandomCentroid(const Pixels& image, int maxIndex, int numberOfCentroids){
    Pixels clusterPoints;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, maxIndex);

    for (int i = 0; i < numberOfCentroids; i++){
        auto pixelIndex = distr(gen);
        clusterPoints.push_back(image[pixelIndex]);
    }
    return clusterPoints;
}

Pixel sumPixelsAssignedToCluster(const std::vector<int>& cluster, const Pixels& image){
    Pixel sum{};
    for (int j = 0;j < cluster.size(); ++j) {
        sum += image[cluster[j]];
    }
    return sum;
}

PixelsSums sumValuesAssignedToCentroids(const std::vector<std::vector<int>>& clusters, const Pixels& image){
    PixelsSums output(clusters.size());;

    for (int i = 0; i < clusters.size(); ++i){
        if (not clusters[i].empty()) {
            output[i].pixel = sumPixelsAssignedToCluster(clusters[i], image);
            output[i].count = clusters[i].size();
        }

    }
    return output;
}

std::vector<std::vector<int>> makeClusters(const Pixels& centroids, const Pixels& image, std::vector<int>& assignment){
    std::vector<std::vector<int>> clusters(centroids.size());

    for (int i = 0; i < image.size(); ++i) {
        auto clusterIndex = assignCluster(image[i], centroids);
        clusters[clusterIndex].push_back(i);
        assignment[i] = clusterIndex;
    }

    return clusters;
}

Pixels reconstructImage(const std::vector<int>& assignment, const Pixels& centroidsPoints){
    Pixels  output;
    for (int i = 0; i < assignment.size(); ++i) {
        int clusterIndex = assignment[i];
        output.push_back(centroidsPoints[clusterIndex]);
    }
    return output;
}

int sendIndexOfDataToSubprocess(int processes, int dataChunkSize, int NumberTag) {
    int start = 0;
    int end = dataChunkSize;
    for (int i = 0; i < processes - 1; ++i){
        MPI_Send(&start, 1, MPI_INT, i + 1, NumberTag, MPI_COMM_WORLD);
        MPI_Send(&end, 1, MPI_INT, i + 1, NumberTag, MPI_COMM_WORLD);
        start += dataChunkSize;
        end += dataChunkSize;
    }
    return start;
}

PixelsSums sumClusterValuesFromSubprocess(PixelsSums sumOfValuesAssignedToCentroids, const PixelsSums& centroidsOfSubprocess) {
    for (int i = 0; i < centroidsOfSubprocess.size(); ++i){
        sumOfValuesAssignedToCentroids[i].pixel += centroidsOfSubprocess[i].pixel;
        sumOfValuesAssignedToCentroids[i].count += centroidsOfSubprocess[i].count;
    }
    return sumOfValuesAssignedToCentroids;
}

PixelsSums sumClusterValuesFromSubprocess(PixelsSums sumOfValuesAssignedToCentroids, int processes, const MpiReceiveParameters& mpiParams) {
     PixelsSums centroids(sumOfValuesAssignedToCentroids.size());
    for (int i = 0; i < processes - 1; ++i){
        MPI_Recv(centroids.data(), centroids.size(), mpiParams.type, i + 1, mpiParams.tag, MPI_COMM_WORLD, mpiParams.status);
        sumOfValuesAssignedToCentroids = sumClusterValuesFromSubprocess(sumOfValuesAssignedToCentroids, centroids);
    }
    return sumOfValuesAssignedToCentroids;
}

Pixels updateCentroidsValues(const PixelsSums& sumOfValuesAssignedToCentroids) {
    Pixels centroidsPoints(sumOfValuesAssignedToCentroids.size());
    for (int i = 0; i < sumOfValuesAssignedToCentroids.size(); ++i){
        if (sumOfValuesAssignedToCentroids[i].count != 0) {
            centroidsPoints[i] = sumOfValuesAssignedToCentroids[i].pixel;
            centroidsPoints[i] /= sumOfValuesAssignedToCentroids[i].count;
        }
    }
    return centroidsPoints;
}

std::vector<int> sumOfAllAssignments(const MpiReceiveParameters& mpiParams, int dataChunkSize, int processes){
    std::vector<int> sumOfAssignments;
    std::vector<int> chunkOfData(dataChunkSize);
    for (int i = 0; i < processes - 1; ++i){
        MPI_Recv(chunkOfData.data(), dataChunkSize, mpiParams.type, i + 1, mpiParams.tag, MPI_COMM_WORLD, mpiParams.status);
        sumOfAssignments.insert(sumOfAssignments.end(), chunkOfData.begin(), chunkOfData.end());
    }
    return sumOfAssignments;
}

std::vector<int> kMeansSubprocessLoop(const KMeansParams& params, const KMeansMpiParams& mpiParams){
    std::vector<int> assignment(params.end-params.start);
    for (int i = 0; i < params.iterations; ++i){
        Pixels centroidsPoints(params.centroids);
        MPI_Bcast(centroidsPoints.data(), params.centroids, mpiParams.pixelDt, 0, MPI_COMM_WORLD);

        Pixels imagePart(params.image.begin() + params.start, params.image.begin() + params.end);
        std::vector<std::vector<int>> clusters = makeClusters(centroidsPoints, imagePart, assignment);
        auto sumOfValuesAssignedToCentroids = sumValuesAssignedToCentroids(clusters, imagePart);
        // MPI_Send(sumOfValuesAssignedToCentroids.data(),
        //          sumOfValuesAssignedToCentroids.size(),
        //          mpiParams.sumPixelDt,
        //          0,
        //          mpiParams.tag,
        //         MPI_COMM_WORLD);
    }
    return assignment;
}

std::pair<std::vector<int>, Pixels> kMeansMainLoop(KMeansParamsFormMainProcess params, const KMeansMpiParams& mpiParams){
    MPI_Status status;
    std::vector<int> assignment(params.image.size());
    for (int i = 0; i < params.iterations; ++i) {
        MPI_Bcast(params.centroidsPoints.data(), params.centroids, mpiParams.pixelDt, 0, MPI_COMM_WORLD);
        std::vector<std::vector<int>> clusters = makeClusters(params.centroidsPoints, params.image, assignment);
        PixelsSums sumOfValuesAssignedToCentroids = sumValuesAssignedToCentroids(clusters, params.image);

        // MpiReceiveParameters paramsForReceive{.type = mpiParams.sumPixelDt, .tag = mpiParams.tag, .status = &status};
        // sumOfValuesAssignedToCentroids = sumClusterValuesFromSubprocess(
        //     sumOfValuesAssignedToCentroids,
        //     params.processes,
        //     paramsForReceive);
        // params.centroidsPoints = updateCentroidsValues(sumOfValuesAssignedToCentroids);
    }
    return {assignment, params.centroidsPoints};
}

int main(int argc, char *argv[]) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        return -1;
    }
    std::string fileName{"img/lennaArray.txt"};
    int width{512};
    int hight{512};
    auto image = getImageFromFile(fileName, width, hight);

    int centroids{std::stoi(argv[1])};
    int iterations{std::stoi(argv[2])};

    int processes;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype pixelDt;
    MPI_Type_contiguous(3, MPI_INT, &pixelDt);
    MPI_Type_commit(&pixelDt);

    MPI_Datatype sumPixelDt;
    MPI_Type_contiguous(2, MPI_INT, &sumPixelDt);
    MPI_Type_commit(&sumPixelDt);
    MPI_Status status;

    const int NumberTag{1};
    const int sumOfPixelsTag{2};
    const int assignmentsTag{3};
    const KMeansMpiParams mpiParams{pixelDt, sumPixelDt, sumOfPixelsTag};
    std::cout << image.size();
    if (rank == 0) {
        const int size = image.size();
        const int dataChunkSize = size / processes;
        Pixels centroidsPoints = getRandomCentroid(image, image.size() - 1, centroids);

        const int start = sendIndexOfDataToSubprocess(processes, dataChunkSize, NumberTag);
        const int end = size;
        const Pixels imagePart(image.begin() + start, image.begin() + end);
        const KMeansParamsFormMainProcess params{centroids, iterations, processes, imagePart, centroidsPoints};

        const auto assignmentsAndCentroidsPoints = kMeansMainLoop(params, mpiParams);
        const std::vector<int> assignment = assignmentsAndCentroidsPoints.first;
        centroidsPoints = assignmentsAndCentroidsPoints.second;


       //const MpiReceiveParameters paramsForReceive{.type = MPI_INT, .tag = assignmentsTag, .status = &status};
        //std::vector<int> sumOfAssignments;// = sumOfAllAssignments(paramsForReceive, dataChunkSize, processes);
        //sumOfAssignments.insert(sumOfAssignments.end(), assignment.begin(), assignment.end());

        imageToFile("img/lennaArray3.txt", reconstructImage(assignment, centroidsPoints));

    } else {
        int start{}, end{};
        MPI_Recv(&start, 1, MPI_INT, 0, NumberTag, MPI_COMM_WORLD, &status);
        MPI_Recv(&end, 1, MPI_INT, 0, NumberTag, MPI_COMM_WORLD, &status);

        KMeansParams params{start, end, centroids, iterations, image};
        std::vector<int> assignment = kMeansSubprocessLoop(params, mpiParams);

        //MPI_Send(assignment.data(), assignment.size(), MPI_INT, 0, assignmentsTag, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}

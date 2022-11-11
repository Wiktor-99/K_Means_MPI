#include <mpi.h>
#include <iostream>
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

std::vector<Pixel> getImageFromFile(std::string fileName, int width, int hight){
    std::fstream file{fileName};
    std::vector<Pixel>  output;
    for (int i = 0; i < width * hight; ++i){
        int r, g, b;
        file >> r >> g >> b;
        output.push_back({r, g, b});
    }
    return output;
}

void imageToFile(std::string fileName, const std::vector<Pixel>& image){
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

int assignCluster(Pixel imagePixel, const std::vector<Pixel>& clusters)
{
    std::vector<int> distances;
	for (const auto& clusterPixel : clusters)
	{
		int distance = euclideanDistance(imagePixel, clusterPixel);
        distances.push_back(distance);

	}
	return std::distance(std::begin(distances), std::min_element(std::begin(distances), std::end(distances)));
}

std::vector<Pixel> getRandomCentroid(const std::vector<Pixel>& image, int maxIndex, int numberOfCentroids){
    std::vector<Pixel> clusterPoints;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, maxIndex);

    for (int i = 0; i < numberOfCentroids; i++){
        auto pixelIndex = distr(gen);
        clusterPoints.push_back(image[pixelIndex]);
    }
    return clusterPoints;
}

Pixel sumPixelsAssignedToCluster(const std::vector<int>& cluster, const std::vector<Pixel>& image){
    Pixel sum{};
    for (int j = 0;j < cluster.size(); ++j) {
        sum += image[cluster[j]];
    }
    return sum;
}

std::vector<PixelSum> updateCentroidValues(const std::vector<std::vector<int>>& clusters, const std::vector<Pixel>& image){
    std::vector<PixelSum> output(clusters.size());;

    for (int i = 0; i < clusters.size(); ++i){
        if (not clusters[i].empty()) {
            output[i].pixel = sumPixelsAssignedToCluster(clusters[i], image);
            output[i].count = clusters[i].size();
        }

    }
    return output;
}

std::vector<std::vector<int>> makeClusters(const std::vector<Pixel>& centroids, const std::vector<Pixel>& image, std::vector<int>& assignment){
    std::vector<std::vector<int>> clusters(centroids.size());

    for (int i = 0; i < image.size(); ++i) {
        auto clusterIndex = assignCluster(image[i], centroids);
        clusters[clusterIndex].push_back(i);
        assignment[i] = clusterIndex;
    }

    return clusters;
}

std::vector<Pixel> reconstructImage(const std::vector<int>& assignment, const std::vector<Pixel>& centroidsPoints){
    std::vector<Pixel>  output;
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

int main(int argc, char *argv[]) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        return -1;
    }
    std::string fileName{"img/lennaArray.txt"};
    int width{512};
    int hight{512};
    auto image = getImageFromFile(fileName, width, hight);
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

    const int NumberTag = 1;
    const int arrayTag = 2;
    const int arrayTag2 = 3;
    MPI_Status status;
    int size = image.size();
    int dataChunkSize = size / processes;

    int centroids{5};
    int iterations{20};

    if (rank == 0) {
        std::vector<Pixel> centroidsPoints = getRandomCentroid(image, image.size() - 1, centroids);

        int start = sendIndexOfDataToSubprocess(processes, dataChunkSize, NumberTag);
        int end = size;
        std::vector<Pixel> imagePart(image.begin() + start, image.begin() + end);
        std::vector<int> assignment(imagePart.size());


        for (int k = 0; k < iterations; ++k) {
            MPI_Bcast(centroidsPoints.data(), centroids, pixelDt, 0, MPI_COMM_WORLD);
            std::vector<std::vector<int>> clusters = makeClusters(centroidsPoints, imagePart, assignment);
            std::vector<PixelSum> centroids2 = updateCentroidValues(clusters, imagePart);

            std::vector<PixelSum> centroids(centroidsPoints.size());
            for (int i = 0; i < processes - 1; ++i){
                MPI_Recv(centroids.data(), centroids.size(), sumPixelDt, i + 1, arrayTag, MPI_COMM_WORLD, &status);
                for (int i = 0; i < centroids.size(); ++i){
                    centroids2[i].pixel += centroids[i].pixel;
                    centroids2[i].count += centroids[i].count;
                }
            }

            for (int i = 0; i < centroids2.size(); ++i){
                centroidsPoints[i] = centroids2[i].pixel;
                centroidsPoints[i] /= centroids2[i].count;
            }

        }


        std::vector<int> chunkOfData(dataChunkSize);
        std::vector<int> sumOfAssignments;
        for (int i = 0; i < processes - 1; ++i){
            MPI_Recv(chunkOfData.data(), dataChunkSize, MPI_INT, i + 1, arrayTag2, MPI_COMM_WORLD, &status);
            sumOfAssignments.insert(sumOfAssignments.end(), chunkOfData.begin(), chunkOfData.end());
        }

        sumOfAssignments.insert(sumOfAssignments.end(), assignment.begin(), assignment.end());

        imageToFile("img/lennaArray3.txt", reconstructImage(sumOfAssignments, centroidsPoints));

    } else {

        int start{}, end{};
        std::vector<int> assignment1;
        MPI_Recv(&start, 1, MPI_INT, 0, NumberTag, MPI_COMM_WORLD, &status);
        MPI_Recv(&end, 1, MPI_INT, 0, NumberTag, MPI_COMM_WORLD, &status);
        for (int i = 0; i < iterations; ++i){
            std::vector<Pixel> centroidsPoints(centroids);
            MPI_Bcast(centroidsPoints.data(), centroids, pixelDt, 0, MPI_COMM_WORLD);

            std::vector<Pixel> imagePart(image.begin() + start, image.begin() + end);
            std::vector<int> assignment(imagePart.size());
            std::vector<std::vector<int>> clusters = makeClusters(centroidsPoints, imagePart, assignment);
            auto centroidsPoints2 = updateCentroidValues(clusters, imagePart);
            MPI_Send(centroidsPoints2.data(), centroidsPoints2.size(), sumPixelDt, 0, arrayTag, MPI_COMM_WORLD);
            assignment1 = assignment;
        }

        MPI_Send(assignment1.data(), assignment1.size(), MPI_INT, 0, arrayTag2, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}

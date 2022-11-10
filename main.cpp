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

int main(){
    std::string fileName{"img/lennaArray.txt"};
    int width{512};
    int hight{512};
    auto image = getImageFromFile(fileName, width, hight);
    int centroids{5};


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, 512*512-1);

    std::vector<Pixel> clusterPoints;


    for (int i = 0; i < centroids; i++){
        auto pixelIndex = distr(gen);
        clusterPoints.push_back(image[pixelIndex]);
    }

    std::vector<int> assignment(image.size());
    for (int iterations = 0; iterations < 500; ++iterations){
        std::vector<std::vector<int>> clusters(centroids);

        for (int i = 0; i < image.size(); ++i) {
            auto clusterIndex = assignCluster(image[i], clusterPoints);
            clusters[clusterIndex].push_back(i);
            assignment[i] = clusterIndex;
        }

        for (int i = 0; i < centroids; ++i){
            for (int j = 0;j < clusters[i].size(); ++j) {
                clusterPoints[i].r += image[clusters[i][j]].r;
                clusterPoints[i].g += image[clusters[i][j]].g;
                clusterPoints[i].b += image[clusters[i][j]].b;
            }
            clusterPoints[i].r /= clusters[i].size();
            clusterPoints[i].g /= clusters[i].size();
            clusterPoints[i].b /= clusters[i].size();
        }
    }

    std::vector<Pixel>  output;

    for (int i = 0; i < width * hight; ++i) {
        int clusterIndex = assignment[i];
        output.push_back(clusterPoints[clusterIndex]);
    }


    imageToFile("img/lennaArray2.txt", output);
}
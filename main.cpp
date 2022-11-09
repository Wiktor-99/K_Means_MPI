#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

struct Image
{
    void appendPixel(int r, int g, int b){
        redPixel.push_back(r);
        greenPixel.push_back(g);
        bluePixel.push_back(b);
    }

    std::vector<int> redPixel;
    std::vector<int> greenPixel;
    std::vector<int> bluePixel;
};

Image getImageFromFile(std::string fileName, int width, int hight){
    std::fstream file{fileName};
    Image output;
    for (int i = 0; i < width * hight; ++i){
        int r, g, b;
        file >> r >> g >> b;
        output.appendPixel(r, g, b);
    }
    return output;
}

void imageToFile(std::string fileName, Image& image){
    std::fstream file{fileName};
    for (int i = 0; i < image.redPixel.size(); ++i){
        file << image.redPixel[i] << ' ' << image.greenPixel[i] << ' ' << image.bluePixel[i] << ' ';
    }
}

double euclideanDistance(int red, int green, int blue, int x, int y, int z){
    return std::sqrt(std::pow(red - x, 2) + std::pow(green - y, 2) + std::pow(blue - z, 2));
}

int assignCluster(int r,int g, int b, const Image& image , int clusters)
{
    std::vector<int> distances;
	for (int i = 0; i < clusters; ++i)
	{
		int distance = euclideanDistance(r, g, b, image.redPixel[i], image.greenPixel[i], image.bluePixel[i]);
        distances.push_back(distance);

	}
	return std::distance(std::begin(distances), std::min_element(std::begin(distances), std::end(distances)));
}

int main(){
    std::string fileName{"img/lennaArray.txt"};
    int width{512};
    int hight{512};
    auto image = getImageFromFile(fileName, width, hight);
    int centroids{7};


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, 512*512-1);

    Image clusterPoints2;


    for (int i = 0; i < centroids; i++){
        auto pixelIndex = distr(gen);
        clusterPoints2.appendPixel(image.redPixel[pixelIndex], image.greenPixel[pixelIndex], image.bluePixel[pixelIndex]);
    }

    std::vector<int> assignment(image.redPixel.size());
    for (int iterations = 0; iterations < 100; ++iterations){
        std::vector<std::vector<int>> clusters(centroids);

        for (int i = 0; i < image.redPixel.size(); ++i) {
            auto clusterIndex = assignCluster(image.redPixel[i], image.greenPixel[i], image.bluePixel[i], clusterPoints2, centroids);
            clusters[clusterIndex].push_back(i);
            assignment[i] = clusterIndex;
        }


        Image clusterPoints;
        clusterPoints.redPixel = std::vector<int>(centroids);
        clusterPoints.greenPixel = std::vector<int>(centroids);
        clusterPoints.bluePixel = std::vector<int>(centroids);

        for (int i = 0; i < centroids; ++i){
            for (int j = 0;j < clusters[i].size(); ++j) {
            clusterPoints.redPixel[i] += image.redPixel[clusters[i][j]];
            clusterPoints.greenPixel[i] += image.greenPixel[clusters[i][j]];
            clusterPoints.bluePixel[i] += image.bluePixel[clusters[i][j]];
            }
            clusterPoints.redPixel[i] /= clusters[i].size();
            clusterPoints.greenPixel[i] /= clusters[i].size();
            clusterPoints.bluePixel[i] /= clusters[i].size();
        }

        clusterPoints2 = clusterPoints;
    }

    Image output;

    for (int i = 0; i < width * hight; ++i) {
        int clusterIndex = assignment[i];
        int red = clusterPoints2.redPixel[clusterIndex];
        int green = clusterPoints2.greenPixel[clusterIndex];
        int blue = clusterPoints2.bluePixel[clusterIndex];

        output.appendPixel(red, green, blue);
    }


    imageToFile("img/lennaArray2.txt", output);
}
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits.h>

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

    imageToFile("img/lennaArray2.txt", image);
}
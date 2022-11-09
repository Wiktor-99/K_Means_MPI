#include <iostream>
#include <fstream>
#include <vector>

struct Image
{
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
        output.redPixel.push_back(r);
        output.greenPixel.push_back(g);
        output.bluePixel.push_back(b);
    }
    return output;
}

void imageToFile(std::string fileName, Image& image){
    std::fstream file{fileName};
    for (int i = 0; i < image.redPixel.size(); ++i){
        file << image.redPixel[i] << ' ' << image.greenPixel[i] << ' ' << image.bluePixel[i] << ' ';
    }
}

int main(){
    std::string fileName{"img/lennaArray.txt"};
    int width{512};
    int hight{512};
    auto image = getImageFromFile(fileName, width, hight);
    imageToFile("img/lennaArray2.txt", image);
}
#include <iostream>
#include <fstream>
#include <vector>

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

int main(){
    std::string fileName{"img/lennaArray.txt"};
    int width{512};
    int hight{512};
    auto image = getImageFromFile(fileName, width, hight);
    imageToFile("img/lennaArray2.txt", image);
}
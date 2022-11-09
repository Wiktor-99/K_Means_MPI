#include <iostream>
#include <fstream>
#include <vector>

struct Pixel
{
    Pixel(int r, int g, int b) : r{r}, g{g}, b{b} {}

    int r;
    int g;
    int b;
};


std::vector<Pixel> getImageFromFile(std::string fileName, int width, int hight){
    std::fstream file{fileName};
    std::vector<Pixel> output;
    for (int i = 0; i < width * hight; ++i){
        int r, g, b;
        file >> r >> g >> b;
        output.emplace_back(r, g, b);
    }
    return output;
}

void imageToFile(std::string fileName, const std::vector<Pixel>& image){
    std::fstream file{fileName};
    for (const auto& [r, g, b] : image){
        file << r << ' ' << g << ' ' << b << ' ';
    }
}

int main(){
    std::string fileName{"img/lennaArray.txt"};
    int width{512};
    int hight{512};
    auto image = getImageFromFile(fileName, width, hight);
    imageToFile("img/lennaArray2.txt", image);
}
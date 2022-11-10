from PIL import Image
import numpy as np
import sys

def imageToFile(pathInput, pathOutput):
    img = Image.open(pathInput)
    array = np.asarray(img, dtype=np.uint8)
    array.tofile(pathOutput, sep=' ')

def fileToArray(inputFile = 'img/lennaArray.txt', width = 512, hight = 512, channels = 3):
    array = np.fromfile(inputFile, dtype=np.uint8, sep= ' ')
    return array[0:width * hight * channels].reshape(width, hight, channels)

def saveArrayAsImage(array, outputName = 'img/lennaAfterClustering.png'):
    image = Image.fromarray(array)
    image.save(outputName, quality=50)


if __name__ == "__main__":
    if len(sys.argv) != 3 and len(sys.argv) != 5:
        print('Wrong number of parameters')
        print('To save image as array, pass path to input file and output file')
        print('To save image from array, pass path to input file, output file, width and hight of image')
    elif len(sys.argv) == 3:
        imageToFile(sys.argv[1], sys.argv[2])
    else:
        array = fileToArray(sys.argv[1], int(sys.argv[3]), int(sys.argv[4]))
        saveArrayAsImage(array, sys.argv[2])
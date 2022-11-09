from PIL import Image
import numpy as np

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
    imageToFile('img/lenna.png', 'img/lennaArray.txt')
    array = fileToArray()
    saveArrayAsImage(array)
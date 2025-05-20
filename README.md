# Basic computer vision algorithms
A lightweight C library implementing core image processing algorithms with CLI support.
## Table of content
  - [Functions](#functions)
  - [Installation](#installation)
  - [Structure](#structure)
  - [Examples](#examples)
  - [References](#references)
  - [Remark](#remark)

## Functions
  - Gaussian blur [-gauss]
  - Median filter [-median]
  - Canny edge detection [-edge]
  - Sharpening [-sharpen]
  - Ridge detection [-ridge]
  - Grayscale conversion [-grayscale]
  - Histogram Equalization [-histogram]
  - Help [-help]
## Installation
Build the project:
   ```bash
   gcc filters.c -shared -fPIC -O2 -o filters.dll
   gcc main.c -o bcva -L. -lfilters -lm
  ```
or simply use c.bat file to compile.

It is nessesary to download `stb_image.h` and `stb_image_write.h` libraries 

## Structure
| File            | Description                                                                |
| ----------------- | ------------------------------------------------------------------ |
| main.c  | CLI interface |
| filters.c | Algorithm implementations |
| filters.h | Header file |
| std/ | mage I/O libraries |
## Examples
1. **Gaussian blur**
   
   ```bash
   ./bcva input.jpg -gauss 5 -o output.jpg
   ```
   ![umbrella (Custom)](https://github.com/user-attachments/assets/3355fa3f-7717-4942-b611-6acf06032a61)
   ![umbrellaout (Custom)](https://github.com/user-attachments/assets/60a2cc17-3fc1-47a1-a64e-a6ed485c577b)

2. **Median filter**

  ```bash
   ./bcva input.jpg -median 3 -o output.jpg
   ```
  ![noise](https://github.com/user-attachments/assets/4390ffb8-fd54-406e-91f6-59ec9cfbd9eb)
  ![noiseout](https://github.com/user-attachments/assets/4a2738f3-67db-479f-9912-9a807e633ef9)

3. **Canny edge detection**

  ```bash
  ./bcva input.jpg -edge 0.5 0.6 -o output.jpg
  ```
  ![car (Custom)](https://github.com/user-attachments/assets/6d37f826-83ae-47df-a3fb-922837b8f122)
  ![carout (Custom)](https://github.com/user-attachments/assets/cdcf541b-fd0b-40fa-ba1c-24ab8b097c46)

4. **Sharpening**

  ```bash
  ./bcva input.jpg -sharpen -o output.jpg
  ```
  ![icon (Custom)](https://github.com/user-attachments/assets/5856daf8-4c6a-49b0-8cce-67caaf099dfc)
  ![iconout (Custom)](https://github.com/user-attachments/assets/19bd6a1e-e2d7-410e-b7ee-22471a836f87)

5. **Ridge detecrion**

  ```bash
  ./bcva input.jpg -ridge -o output.jpg
  ```
  ![sandwich (Custom)](https://github.com/user-attachments/assets/d190de99-4fa5-4392-af62-d0f00337593f)
  ![sandwichout (Custom)](https://github.com/user-attachments/assets/39927195-216a-4894-97e1-49c52f56edc9)

6. **Grayscale conversion**

  ```bash
  ./bcva input.jpg -ridge -o output.jpg
  ```
  ![clown (Custom)](https://github.com/user-attachments/assets/be55f631-b1cc-4097-a695-5223cee3a1df)
  ![clownout (Custom)](https://github.com/user-attachments/assets/676052f4-d9de-407b-8cbd-e9fa9a1ac2f4)

7. **Histogram Equalization**

  ```bash
  ./bcva input.jpg -histogram -o output.jpg
  ```
![cat](https://github.com/user-attachments/assets/4919a834-cd48-4a09-ade5-8becbc6cef2d)
![catout](https://github.com/user-attachments/assets/234fb3e9-eaeb-432b-8796-b739da7bdbdf)

## References
  - [opencv](https://github.com/opencv/opencv)
  - [Gaussian blur on wikipedia](https://ru.wikipedia.org/wiki/%D0%A0%D0%B0%D0%B7%D0%BC%D1%8B%D1%82%D0%B8%D0%B5_%D0%BF%D0%BE_%D0%93%D0%B0%D1%83%D1%81%D1%81%D1%83)
  - [Быстрое размытие по Гауссу by vladimirovich on habr](https://habr.com/ru/articles/151157/)
  - [Median filter on wikipedia](https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D0%B4%D0%B8%D0%B0%D0%BD%D0%BD%D1%8B%D0%B9_%D1%84%D0%B8%D0%BB%D1%8C%D1%82%D1%80)
  - [Детектор границ Канни by impersonalis on habr](https://habr.com/ru/articles/114589/)
  - [Комментарий к алгоритму выделения контуров Канни by AndreyIvanoff on habr](https://habr.com/ru/articles/114766/)
  - [Convolution on wikipedia](https://ru.wikipedia.org/wiki/%D0%AF%D0%B4%D1%80%D0%BE_%D1%81%D0%B2%D0%B5%D1%80%D1%82%D0%BA%D0%B8)
  - [F#, увеличение контраста изображения путем выравнивания гистограмм by zverjuga on habr](https://habr.com/ru/articles/502986/)
  - [stb_image and stb_image_write libraries source](https://github.com/nothings/stb)

## Remark

Novosibirsk State University (NSU)
Bachelor's Program: 15.03.06 - Mechatronics and Robotics (AI Profile)

*Imperative Programming Course - Graduation Project*

*Computer vision algorithms in C*




   


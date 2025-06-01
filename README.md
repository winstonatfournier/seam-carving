# seam-carving

This project implements a seam carving algorithm for content-aware image resizing in C.


## Features

- Grayscale energy map calculation
- Dynamic programming approach for seam cost evaluation
- Recursive backtracking to compute all vertical seams
- Image data handling using raw raster access
- Raw memory management via quadruple pointers


## Notes

- `c_img.c` was provided (for reading/writing PPM images and defining `struct rgb_img`). All image processing and seam carving logic was otherwise implemented by myself.
- All helper functions individually passed the black-box autograder. Written under tight deadlines; performance and cleanliness could be improved.


> Development period: April 2024
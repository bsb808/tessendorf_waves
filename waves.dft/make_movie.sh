#!/bin/bash
#avconv -framerate 4 -i "./movie/video%04d.tga" -c:v libx264 -r 30 -pix_fmt yuv420p -crf 20 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" out.mp4


avconv -framerate 30 -i "./movie/video%04d.tga" -c:v h264 -crf 1 out.mov

#!/usr/bin/env python3

# From https://learn.astropy.org/tutorials/1-SpectroscopicTraceTutorial.html

from PIL import Image
import numpy as np

import matplotlib.pyplot as plt

# Step 1: Examine the spectrum

plt.rcParams["image.origin"] = "lower"
plt.style.use("dark_background") # Optional!

spectrum_filename = "aldebaran_3s_1.bmp"

image_data = Image.open(spectrum_filename)

# our data are unsigned 8-bit integers (0-255) representing a monochromatic image
# we can see this by printing the array version of the image
# we can also see its shape, verifying that it is indeed 2-dimensional
image_array = np.array(image_data)
print(image_array)
print(image_array.shape)

# but we"d like to see it with axes labeled
plt.imshow(image_data, cmap="gray")
plt.colorbar(); # the semicolon at the end of the last line prevents ipython from printing out the object

plt.savefig("aldebaran_3s_1.png")

# Step 2a. Try to find the spine to trace using argmax

yvals = np.argmax(image_data, axis=0)
xvals = np.arange(image_data.width)
plt.plot(xvals, yvals, 'x')
plt.ylabel("Argmax trace data")
plt.xlabel("X position");
plt.savefig("spine_mask.png")

bad_pixels = (yvals < 400) | (yvals > 500)

plt.plot(xvals, yvals, 'x')
plt.plot(xvals[bad_pixels], yvals[bad_pixels], 'rx')
plt.ylabel("Argmax trace data")
plt.xlabel("X position");
plt.savefig("spine_mask_bad_pixels.png")

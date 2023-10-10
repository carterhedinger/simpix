# simpix
This program takes the pixels of one image and rearranges them to look like a 2nd image.
This is done through a clever use of the process of annealing in Statistical Physics by
modeling the difference in the color of the pixels in each image as a sort of Ising model.
For accurate results, this program needs to run for several hours, because annealing is a
slow process and typically the slower it is the more accurate results are achieved. For this
reason, my final reasults were calculated using slurm on the Rivanna supercluster.
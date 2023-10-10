======================================== README ========================================

>> Description 

    This program takes the pixels of one image and rearranges them to look like a 2nd image.
    This is done through a clever use of the process of annealing in Statistical Physics by
    modeling the difference in the color of the pixels in each image as a sort of Ising model.
    For accurate results, this program needs to run for several hours, because annealing is a
    slow process and typically the slower it is the more accurate results are achieved. For this
    reason, my final reasults were calculated using slurm on the Rivanna supercluster.

>> Compiling and Running

    This program uses a Makefile to compile. The program should be ran using the following syntax:
    ./simpix SourceFileName.png TargetFileName.png OutputFileName.png
    If you have access to a supercluster, then I would suggest running it using a slurm file to run it.

>> Output

    Both final_out1.png and final_out2.png are both the final results of running this program on the
    two input images Claude_Monet_Le_Grand_Canal.png and Claude_Monet_The_Cliffs_at_Etretat.png, each
    with one as the "source" image and the other as the "target" image. Both final_collage1.png and
    final_collage2.png contain three images each. The top left is the "source" image, the top right is
    the "target" image, and the bottom left is the output image i.e. the "source" image's pixels 
    rearranged to look like the "target" image.
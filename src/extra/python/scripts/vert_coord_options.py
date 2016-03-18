import numpy as np
import matplotlib.pyplot as plt 


def even_sigma_calc(num_levels):
	"The even sigma calculation just divides the atmosphere up into equal sigma increments between 1 and 0. So the height of the model is really set by your number of levels, as the higher the number of levels you have, the smaller your increment will be between 0hPa at the top and whatever your next level down is."
	num_levels_calc = num_levels -1

	b=np.zeros(num_levels)

	for k in np.arange(1,num_levels):
		b[k-1] = (k-1.)/num_levels_calc
	
	b[num_levels-1]=1.

	p_half = b

	flipped_p_half = p_half[::-1] #Note that this is the opposite convention to the model, but done for visual clarity in plots. 

	return flipped_p_half


def uneven_sigma_calc(num_levels, surf_res, exponent, scale_heights):
	"The uneven sigma calculation first splits up the atmosphere into equal increments between 0 and 1, and then does different vertical spacings depending on the parameters. For example, if surf_res = 1 then you get even spacing in height. If surf_res = 0 then you get a height depending on zeta**exponent. For surf_res in between, you get a mix of the two. scale_heights sets the model top height, and exponent determines how heavily biased your level spacings are towards the troposphere. Larger exponent values are more tropospherically biased. "
	num_levels_calc = num_levels -1

	b=np.zeros(num_levels)

	for k in np.arange(1,num_levels):
		zeta = 1. - (k-1.)/num_levels_calc

		Z = surf_res*zeta + (1. - surf_res)*np.power(zeta,exponent)

		b[k-1] = np.exp(-Z*scale_heights)
	
	b[num_levels-1]=1.
#	b[0]=0.


	p_half = b

	flipped_p_half = p_half[::-1] #Note that this is the opposite convention to the model, but done for visual clarity in plots. 

	return flipped_p_half


def p_half_to_p_full(p_half, num_levels):

	p_full = np.zeros(num_levels-1)

	# 0 to num_levels-1 is so that we go through all p_half[k], but we can't have p_half[k+1] with a top number of num_levels-1, so must be num_levels-2. BUT because np.arange doesn't include the end point, we use num_levels-1. 
	for k in np.arange(0,num_levels-1):
		alpha  = 1.0  - p_half[k]*( np.log(p_half[k+1]) - np.log(p_half[k]) )/ (p_half[k+1] - p_half[k])
		p_full[k] = p_half[k+1] * np.exp(-alpha)

	return p_full


if __name__ == "__main__":

	num_levels = 41 #number of half levels

	vert_coord_option = "uneven_sigma"
	surf_res                = 0.5
	scale_heights           = 11.0
	exponent                = 7.0

	if vert_coord_option == "uneven_sigma":

		plt.figure(1)
		for surf_res in [0., 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]:
		
			p_half = uneven_sigma_calc(num_levels, surf_res, exponent, scale_heights)
			p_full = p_half_to_p_full(p_half, num_levels)
			plt.plot(-np.log(p_half), label=surf_res)

		plt.legend(loc='upper left')
		plt.title('Uneven sigma, varying surf_res')

		plt.figure(2)

		for scale_heights in [6., 8., 10., 12.]:
			
			surf_res = 1.0

			p_half = uneven_sigma_calc(num_levels, surf_res, exponent, scale_heights)
			p_full = p_half_to_p_full(p_half, num_levels)
			plt.plot(-np.log(p_half), label=scale_heights)

		plt.legend(loc='upper left')
		plt.title('Uneven sigma, varying scale_heights')

		plt.figure(3)

		for exponent in [2., 4., 6., 7., 8., 10., 12.]:
			
			surf_res = 0.
			scale_heights = 11.

			p_half = uneven_sigma_calc(num_levels, surf_res, exponent, scale_heights)
			p_full = p_half_to_p_full(p_half, num_levels)
			plt.plot(-np.log(p_half), label=exponent)

		plt.legend(loc='upper left')
		plt.title('Uneven sigma, varying exponent')


		plt.show()


	if vert_coord_option == "even_sigma":

		p_half = even_sigma_calc(num_levels)
		p_full = p_half_to_p_full(p_half, num_levels)
		plt.plot((p_half))
		plt.show()





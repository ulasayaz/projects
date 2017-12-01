from __future__ import print_function
import argparse
import os
import random
import torch
import torch.nn as nn
import torch.nn.parallel
import torch.backends.cudnn as cudnn
import torch.optim as optim
import torch.utils.data
import torchvision.datasets as dset
import torchvision.transforms as transforms
import torchvision.utils as vutils
from torch.autograd import Variable

import matplotlib.pyplot as plt
import numpy as np
import time

from scipy.optimize import minimize, rosen, rosen_der
# Refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
from models import _netG
import dataloader as nyu_dset

parser = argparse.ArgumentParser(description='Image Reconstruction')
parser.add_argument('--dataset', default='nyudepthv2', help='nyudepthv2 | fake')
parser.add_argument('--cuda', action='store_true', help='enables cuda')
parser.add_argument('--modality', default='rgb', help='rgb | rgbd')
parser.add_argument('-s', '--sampling-strategy', default='rg', help="options: uniform, or any combination of 'r', 'g' and 'b'")
parser.add_argument('-p', '--percentage-samples', type=float, default=1, help='percentage of samples')
opt = parser.parse_args()

if torch.cuda.is_available() and not opt.cuda:
	print("WARNING: You have a CUDA device, so you should probably run with --cuda")

image_size = 128
nz = 100 # size of the latent z vector
ngf = 64
nc = 3
N = nc * image_size * image_size

criterion = nn.MSELoss(size_average=False)
# criterion = nn.L1Loss(size_average=False)

if opt.sampling_strategy == 'uniform':
	sampling_matrix = (np.random.uniform(0.0, 1.0, size=(1,image_size,image_size)) < opt.percentage_samples).astype(float)
	sampling_matrix = np.repeat(sampling_matrix, nc, axis=0)
else:
	opt.percentage_samples = 0
	if 'r' in opt.sampling_strategy:
		sampling_matrix_r = np.ones(shape=(1,image_size,image_size))
		opt.percentage_samples += 1.0/3
	else:
		sampling_matrix_r = np.zeros(shape=(1,image_size,image_size))

	if 'g' in opt.sampling_strategy:
		sampling_matrix_g = np.ones(shape=(1,image_size,image_size))
		opt.percentage_samples += 1.0/3
	else:
		sampling_matrix_g = np.zeros(shape=(1,image_size,image_size))

	if 'b' in opt.sampling_strategy:
		sampling_matrix_b = np.ones(shape=(1,image_size,image_size))
		opt.percentage_samples += 1.0/3
	else:
		sampling_matrix_b = np.zeros(shape=(1,image_size,image_size))
	sampling_matrix = np.concatenate((sampling_matrix_r,sampling_matrix_g,sampling_matrix_b), axis=0)

sampling_matrix_pth = Variable(torch.from_numpy(sampling_matrix).float(), volatile=False)
identity_matrix_pth = Variable(torch.ones(nc, image_size, image_size).float(), volatile=False)
if opt.cuda:
	sampling_matrix_pth = sampling_matrix_pth.cuda()
	identity_matrix_pth = identity_matrix_pth.cuda()
# print(sampling_matrix.shape)
# print(float(np.sum(sampling_matrix))/sampling_matrix.size)

iter = 0
def empirical_risk(netG, x, y, A):
	# x is numpy vector, y is torch variable
	global iter
	iter = iter+1

	# convert x from a numpy vector to a 4-dim torch float tensor
	x = torch.from_numpy(x).float().view(1, nz, 1, 1)
	if opt.cuda:
		x = x.cuda()

	v = Variable(x, requires_grad = True)
	output = netG(v) * A

	# compute the reconstruction error
	loss = criterion(output, y)
	loss_np = loss.data.cpu().numpy()[0]

	if iter % 50 == 0:
		print('==> iter={}, loss={}'.format( iter, loss_np))

	return loss_np

def reconstruction_error(netG, x, y):
	# x is numpy vector, y is torch variable

	# convert x from a numpy vector to a 4-dim torch float tensor
	x = torch.from_numpy(x).float().view(1, nz, 1, 1)
	if opt.cuda:
		x = x.cuda()

	v = Variable(x, requires_grad = False)
	output = netG(v)

	# compute the reconstruction error
	loss = criterion(output, y)
	loss_np = loss.data.cpu().numpy()[0]

	return loss_np

def grad(netG, x, y, A):
	# convert x from a numpy vector to a 4-dim torch float tensor
	x = torch.from_numpy(x).float().view(1, nz, 1, 1)
	if opt.cuda:
		x = x.cuda()

	v = Variable(x, requires_grad = True)
	output = netG(v) * A

	# compute the reconstruction error
	loss = criterion(output, y)
	loss_np = loss.data.cpu().numpy()[0]

	# compute gradient
	loss.backward()
	grad_np = np.squeeze(v.grad.data.cpu().numpy())

	return grad_np

# def optimize(netG, x0, y):
# 	# TODO: gradient descent on GPU

# 	x = torch.from_numpy(x0).float().view(1, nz, 1, 1)
# 	if opt.cuda:
# 		x = x.cuda()
# 	v = Variable(x, requires_grad = True)

# 	for i in range(10000):
# 		output = netG(v)

# 		# compute the reconstruction error
# 		loss = criterion(output, y)
# 		loss_np = loss.data.cpu().numpy()[0]

# 		# compute gradient
# 		loss.backward()

# 		# update x
# 		step_size = 1e4
# 		v.data -= step_size * v.grad.data
# 		v.grad.data.zero_()
# 		print('iter={}, loss={}'.format(i, loss.data.cpu().numpy()[0]))

# 	return v.data

def main():
	global iter

	# cuda
	if opt.cuda:
		# torch.cuda.manual_seed_all(opt.manualSeed)
		netG = _netG(1)
		netG.cuda()
		criterion.cuda()
	else:
		netG = _netG(0)
	netG.load_state_dict(torch.load('results/nyudepthv2.modality=rgb.nz=100.ngf=64.bs=64/netG_epoch_55.pth'))
	print(netG)

	if opt.dataset == 'nyudepthv2':
		# create dataloader for real image test
		valdir = os.path.join('data', 'nyudepthv2', 'val')
		dataset = nyu_dset.ImageFolder(valdir, type='val', modality=opt.modality)
		input_tensor, depth_tensor = dataset.__getitem__(10)
		if opt.cuda:
			input_tensor = input_tensor.cuda()
			depth_tensor = depth_tensor.cuda()
		y_gt = Variable(input_tensor, volatile=False)
		y = sampling_matrix_pth * y_gt
	elif opt.dataset == 'fake':
		# create random ground truth signal of interest
		x_gt_np = np.random.randn(nz)
		x_gt_pth = torch.from_numpy(x_gt_np).float().view(1, nz, 1, 1)
		if opt.cuda:
			x_gt_pth = x_gt_pth.cuda()
		y_gt = netG(Variable(x_gt_pth, volatile=True))
		y = sampling_matrix_pth * netG(Variable(x_gt_pth, volatile=False))
	
	y = y.detach()

	start_time = time.time()

	# optimization with gradient
	print('==> running solver..')
	iter = 0
	x0_np = np.random.randn(nz)
	result = minimize(lambda x: empirical_risk(netG=netG, x=x, y=y, A=sampling_matrix_pth), 
		x0_np, method='BFGS', jac=lambda x: grad(netG=netG, x=x, y=y, A=sampling_matrix_pth))

	# evaluate reconstruction error without subsampling
	l2norm_gt = torch.sum(y_gt.data * y_gt.data)
	perc_error = reconstruction_error(netG=netG, x=result.x, y=y_gt) / l2norm_gt

	# customized gradient descent algorithm
	# res = optimize(netG, x0_np, y)
	
	print('==> pixel-wise average reconstruction error={:.2f}%'.format(100*perc_error))
	print('==> reconstruction time={:.2f}s'.format(time.time()-start_time))

	# visualize reconstruction and ground truth
	# Two subplots, the axes array is 1-d
	fig = plt.figure()
	ax1 = fig.add_subplot(1,3,1, adjustable='box', aspect=1)
	ax2 = fig.add_subplot(1,3,2, adjustable='box', aspect=1)
	ax3 = fig.add_subplot(1,3,3, adjustable='box', aspect=1)

	# measurements
	img_measurements = np.squeeze(y.data.cpu().numpy()) / 2 + 0.5
	img_measurements = np.transpose(img_measurements * sampling_matrix, (1, 2, 0)) 
	ax1.imshow(img_measurements, extent=[0,100,0,1], aspect=100)
	ax1.set_title('samples={:.2f}%'.format(100*opt.percentage_samples))
	ax1.axis('off')

	# ground truth
	img_gt = np.transpose(np.squeeze(y_gt.data.cpu().numpy()), (1, 2, 0))
	ax2.imshow(img_gt/2 + 0.5, extent=[0,100,0,1], aspect=100)
	ax2.set_title('ground truth')
	ax2.axis('off')

	if opt.cuda:
	    v_rec = Variable(torch.from_numpy(result.x).cuda().float().view(1, nz, 1, 1))
	    rec_pth = netG(v_rec).data
	    img_rec = np.transpose(np.squeeze(rec_pth.cpu().numpy()), (1, 2, 0))
	else:
	    v_rec = Variable(torch.from_numpy(result.x).float().view(1, nz, 1, 1))
	    rec_pth = netG(v_rec).data
	    img_rec = np.transpose(np.squeeze(rec_pth.numpy()), (1, 2, 0))
	ax3.imshow(img_rec/2 + 0.5, extent=[0,100,0,1], aspect=100)
	ax3.set_title( 'error={:.2f}%'.format(100*perc_error) )
	ax3.axis('off')

	# save results
	fig.savefig('reconstruction.pdf', format='pdf')

	plt.show()

if __name__ == '__main__':
	main()


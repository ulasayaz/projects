import argparse
import torch
import torch.nn as nn
import torch.backends.cudnn as cudnn
cudnn.benchmark = True
import torch.optim as optim
from models import *

def weights_init(m):
    # Initialize filters with Gaussian random weights
    if isinstance(m, nn.Conv2d):
        n = m.kernel_size[0] * m.kernel_size[1] * m.out_channels
        m.weight.data.normal_(0, math.sqrt(2. / n))
        if m.bias is not None:
            m.bias.data.zero_()
    elif isinstance(m, nn.ConvTranspose2d):
        n = m.kernel_size[0] * m.kernel_size[1] * m.in_channels
        m.weight.data.normal_(0, math.sqrt(2. / n))
        if m.bias is not None:
            m.bias.data.zero_()
    elif isinstance(m, nn.BatchNorm2d):
        m.weight.data.fill_(1)
        m.bias.data.zero_()

assert torch.cuda.is_available(), "Error: CUDA not available"

parser = argparse.ArgumentParser(description='Invertible Neural Networks')
parser.add_argument('-i', '--iterations', default=1000, help='number of iterations for backprop')
parser.add_argument('--lr', '--learning-rate', default=0.01, type=float, 
    metavar='LR', help='learning rate (default 0.01)')
parser.add_argument('-ld', '--latent-dim', default=100, help='dimension of latent vector')
args = parser.parse_args()

criterion = nn.MSELoss(reduction='sum')
model = _netG().cuda()
model.apply(weights_init) # Gaussian random weights
model.eval()

z_gt = torch.randn(1, args.latent_dim, 1, 1).cuda()  # ground truth target latent code
y_gt = model(z_gt) # the corresponding ground truth image
print("latent size = {}, output image size = {}".format(z_gt.size(), y_gt.size()))

z_estimate = torch.randn(1, args.latent_dim, 1, 1).cuda() # our estimate, initialized randomly
z_estimate.requires_grad = True
optimizer = optim.Adam([z_estimate], lr=args.lr)
for i in range(args.iterations):
    y_estimate = model(z_estimate)
    optimizer.zero_grad()
    loss = criterion(y_estimate, y_gt.detach())
    if i % 20 == 0:
        print("iter {:04d}: y_error = {:03g}, z_error={:03g}".format(i,
            loss.item(), criterion(z_estimate, z_gt.detach())))
    loss.backward()
    optimizer.step()

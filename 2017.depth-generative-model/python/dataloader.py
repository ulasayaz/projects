# modified from https://github.com/pytorch/vision/blob/master/torchvision/datasets/folder.py


import os
import os.path
import numpy as np
import torch.utils.data as data
import h5py
import transforms

IMG_EXTENSIONS = [
    '.h5',
]

def is_image_file(filename):
    return any(filename.endswith(extension) for extension in IMG_EXTENSIONS)


def find_classes(dir):
    classes = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    classes.sort()
    class_to_idx = {classes[i]: i for i in range(len(classes))}
    return classes, class_to_idx


def make_dataset(dir, class_to_idx):
    images = []
    dir = os.path.expanduser(dir)
    for target in sorted(os.listdir(dir)):
        # print(target)
        d = os.path.join(dir, target)
        if not os.path.isdir(d):
            continue

        for root, _, fnames in sorted(os.walk(d)):
            for fname in sorted(fnames):
                if is_image_file(fname):
                    path = os.path.join(root, fname)
                    item = (path, class_to_idx[target])
                    images.append(item)

    return images

def h5_loader(path):
    h5f = h5py.File(path, "r")
    rgb = np.array(h5f['rgb'])
    rgb = np.transpose(rgb, (1, 2, 0))
    depth = np.array(h5f['depth'])

    return rgb, depth


iheight, iwidth = 480, 640
# oheight, owidth = 228, 304
image_size = 128
oheight, owidth = image_size, image_size
color_jitter = transforms.ColorJitter(0.1, 0.1, 0.1)
# ImageNet
# normalize = transforms.Normalize(mean = [ 0.485, 0.456, 0.406 ], std = [ 0.229, 0.224, 0.225 ])
normalize = transforms.NormalizeNumpyArray(mean = [ 0.5, 0.5, 0.5 ], std = [ 0.5, 0.5, 0.5 ])

def train_transform(rgb, depth):
    s = np.random.uniform(1.0, 1.5) # random scaling
    # print("scale factor s={}".format(s))
    depth_np = depth / s
    angle = np.random.uniform(-5.0, 5.0) # random rotation degrees
    do_flip = np.random.uniform(0.0, 1.0) < 0.5 # random horizontal flip

    # set zeros in depth as NaN
    depth_np[depth_np == 0] = np.nan

    # perform 1st part of data augmentation
    transform = transforms.Compose([
        transforms.Resize(float(image_size) / iheight), # this is for computational efficiency, since rotation is very slow
        transforms.Rotate(angle),
        transforms.Resize(s),
        transforms.CenterCrop((oheight, owidth)),
        transforms.HorizontalFlip(do_flip),
    ])
    rgb_np = transform(rgb)

    # random color jittering 
    rgb_np = color_jitter(rgb_np)

    rgb_np = np.asfarray(rgb_np, dtype='float') / 255
    rgb_np = normalize(rgb_np) # from [0,1] to [-1,1]

    depth_np = transform(depth_np)
    depth_np[np.isnan(depth_np)] = 0

    depth_np = depth_np / 10.0

    return rgb_np, depth_np

def val_transform(rgb, depth):
    # perform 1st part of data augmentation
    transform = transforms.Compose([
        transforms.Resize(float(image_size) / iheight),
        transforms.CenterCrop((oheight, owidth)),
    ])
    rgb_np = transform(rgb)
    rgb_np = np.asfarray(rgb_np, dtype='float') / 255
    rgb_np = normalize(rgb_np) # from [0,1] to [-1,1]

    depth_np = transform(depth)
    depth_np = depth_np / 10.0

    return rgb_np, depth_np

def rgb2grayscale(rgb):
    return rgb[:,:,0] * 0.2989 + rgb[:,:,1] * 0.587 + rgb[:,:,2] * 0.114


# NYU Depth V2
# nyu_mean = [ 0.48762784,0.41758214,0.40050098 ]
# nyu_std = [ 0.26801008,0.27798527,0.28912789 ]
# normalize_tensor = transforms.NormalizeTensor(mean = nyu_mean, std = nyu_std)
# normalize_np = transforms.NormalizeNumpyArray(mean = nyu_mean, std = nyu_std)
to_tensor = transforms.ToTensor()

class ImageFolder(data.Dataset):
    """A generic data loader where the images are arranged in this way: 

        root/dog/xxx.h5
        root/dog/xxy.h5
        root/dog/xxz.h5

        root/cat/123.h5
        root/cat/nsdf3.h5
        root/cat/asd932_.h5

    Args:
        root (string): Root directory path.
        transform (callable, optional): A function/transform that  takes in both RGB and Depth
            and returns a transformed version. E.g, ``transforms.RandomCrop``
        loader (callable, optional): A function to load an image given its path.

     Attributes:
        classes (list): List of the class names.
        class_to_idx (dict): Dict with items (class_name, class_index).
        imgs (list): List of (image path, class_index) tuples
    """
    modality_names = ['rgb', 'rgbd'] 

    def __init__(self, root, type, modality='rgb', loader=h5_loader):
        classes, class_to_idx = find_classes(root)
        imgs = make_dataset(root, class_to_idx)
        if len(imgs) == 0:
            raise(RuntimeError("Found 0 images in subfolders of: " + root + "\n"
                               "Supported image extensions are: " + ",".join(IMG_EXTENSIONS)))

        self.root = root
        self.imgs = imgs
        self.classes = classes
        self.class_to_idx = class_to_idx
        if type == 'train':
            self.transform = train_transform
        elif type == 'val':
            self.transform = val_transform
        else:
            raise (RuntimeError("Invalid dataset type: " + type + "\n"
                                "Supported dataset types are: train, val"))
        self.loader = loader

        if modality in self.modality_names:
            self.modality = modality
        else:
            raise (RuntimeError("Invalid modality type: " + modality + "\n"
                                "Supported dataset types are: " + ''.join(self.modality_names)))

    def create_rgbd(self, rgb, depth):
        rgbd = np.append(rgb, np.expand_dims(depth, axis=2), axis=2)
        return rgbd

    def __getraw__(self, index):
        """
        Args:
            index (int): Index

        Returns:
            tuple: (rgb, depth) the raw data.
        """
        path, target = self.imgs[index]
        rgb, depth = self.loader(path)
        return rgb, depth

    def __get_all_item__(self, index):
        """
        Args:
            index (int): Index

        Returns:
            tuple: (image, target) where target is the depth image.
        """
        rgb, depth = self.__getraw__(index)
        if self.transform is not None:
            rgb_np, depth_np = self.transform(rgb, depth)
        else:
            raise(RuntimeError("transform not defined"))

        # color normalization
        # rgb_tensor = normalize_rgb(rgb_tensor)
        # rgb_np = normalize_np(rgb_np)
        
        if self.modality == 'rgb':
            input_np = rgb_np
        elif self.modality == 'rgbd':
            input_np = self.create_rgbd(rgb_np, depth_np)

        input_tensor = to_tensor(input_np)
        depth_tensor = to_tensor(depth_np)
        depth_tensor = depth_tensor.unsqueeze(0)

        return input_tensor, depth_tensor, input_np, depth_np

    def __getitem__(self, index):
        """
        Args:
            index (int): Index

        Returns:
            tuple: (image, target) where target is the depth image.
        """
        input_tensor, depth_tensor, input_np, depth_np = self.__get_all_item__(index)

        return input_tensor, depth_tensor

    def __len__(self):
        return len(self.imgs)
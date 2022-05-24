# invert-convolutional-generative-networks
This repo implements the inversion a convolutional generative neural network for the NeurIPS (formerly abbreviated as NIPS) 2018 paper ["Invertibility of Convolutional Generative Networks from Partial Measurements"](https://papers.nips.cc/paper/2018/hash/e0ae4561193dbf6e4cf7e8f4006948e3-Abstract.html) by Fangchang Ma, Ulas Ayaz, and Sertac Karaman at MIT.

## dcgan
PyTorch code for training a Generative Adversarial Network (GAN) model

## matlab
MatLab code for testing image reconstruction based on a generative model (neural network).
- The Matlab neural network toolbox seems very limited. It is impossible to construct a generative model which takes a vector input.

## python
Python code for testing image reconstruction based on a generative model (neural network).

#### prerequisite
- install [PyTorch](http://pytorch.org/)
- install the python scipy package with `pip install scipy`

#### usage
- If you have a Nvidia GPU and Cuda installed, run `python main.py --cuda`
- Otherwise, run on cpu with `python main.py` 


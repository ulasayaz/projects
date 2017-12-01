# 2018.generative-model
iterative depth reconstruction using generative model and gradient descent

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
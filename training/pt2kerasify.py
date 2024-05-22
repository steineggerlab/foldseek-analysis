#!/usr/bin/env python3

import sys

import torch
import keras

import kerasify

input_fn = sys.argv[1]
output_fn = sys.argv[2]

# Load PyTorch model
encoder = torch.load(input_fn)

# Convert into keras model
# - only ReLU
# - only fully connected
# - no batch norm etc.
input_shape = [encoder[0].in_features]
n_units_lst = [layer.out_features for layer in encoder if not isinstance(layer, torch.nn.ReLU)]

model = keras.models.Sequential()
model.add(keras.layers.Dense(n_units_lst[0], input_shape=input_shape, activation='relu'))
for n_units in n_units_lst[1:-1]:
    model.add(keras.layers.Dense(n_units, activation='relu'))
model.add(keras.layers.Dense(n_units_lst[-1], activation='linear'))

print()
model.summary()

# Copy weights
n_layers = len(n_units_lst)
for i in range(n_layers):
    model.layers[i].set_weights([encoder[i * 2].weight.detach().T, encoder[i * 2].bias.detach()])

# Kerasify
kerasify.export_model(model, output_fn)

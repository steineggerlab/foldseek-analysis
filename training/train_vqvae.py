#! /usr/bin/env python3

import sys
import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

class Decoder(nn.Module):
    def __init__(self, input_dim, hidden_dim, z_dim):
        super(Decoder, self).__init__()
        self.layers = nn.Sequential(nn.Linear(z_dim, hidden_dim),
                        nn.BatchNorm1d(hidden_dim), nn.ReLU(),
                        nn.Linear(hidden_dim, hidden_dim),
                        #nn.BatchNorm1d(hidden_dim), nn.ReLU()
                        )

        self.mu = nn.Linear(hidden_dim, input_dim)
        self.logvar = nn.Linear(hidden_dim, input_dim)

    def forward(self, x):
        x = self.layers(x)
        var = torch.exp(0.5 * self.logvar(x))  # ...should delete 0.5 (NN learns factor anyway)
        return self.mu(x), var


class VectorQuantizer(nn.Module):
    def __init__(self, N_STATES, Z_DIM):
        super(VectorQuantizer, self).__init__()

        self.embedding = nn.Embedding(N_STATES, Z_DIM)
        self.embedding.weight.data.uniform_(-1/N_STATES, 1/N_STATES)

        self.COMMITMENT_COST = 0.25
        self.N_STATES = N_STATES

    def forward(self, inputs):
        # Calculate distances
        distances = (torch.sum(inputs**2, dim=1, keepdim=True)
                    + torch.sum(self.embedding.weight**2, dim=1)
                    - 2 * torch.matmul(inputs, self.embedding.weight.t()))

        # Encoding
        encoding_indices = torch.argmin(distances, dim=1).unsqueeze(1)
        encodings = torch.zeros(encoding_indices.shape[0], self.N_STATES, device=inputs.device)
        encodings.scatter_(1, encoding_indices, 1)

        # Quantize and unflatten
        quantized = torch.matmul(encodings, self.embedding.weight)

        # Loss
        e_latent_loss = F.mse_loss(quantized.detach(), inputs)
        q_latent_loss = F.mse_loss(quantized, inputs.detach())
        loss = q_latent_loss + self.COMMITMENT_COST * e_latent_loss

        quantized = inputs + (quantized - inputs).detach()
        avg_probs = torch.mean(encodings, dim=0)
        perplexity = torch.exp(-torch.sum(avg_probs * torch.log(avg_probs + 1e-10)))

        return loss, quantized, perplexity, encodings


class VAE_VQ(nn.Module):
    def __init__(self, encoder, decoder, z_dim, n_states):
        super(VAE_VQ, self).__init__()

        self.encoder = encoder
        self.decoder = decoder
        self.vq = VectorQuantizer(n_states, z_dim)

    def forward(self, x):
        z = self.encoder(x)
        loss, quantized, perplexity, encodings = self.vq(z)
        mu, var = self.decoder(quantized)
        return loss, (mu, var), perplexity, encodings


def batched(data, batch_size=False):
    if batch_size:
        return [[elem[i * batch_size:(i + 1) * batch_size] for elem in data]
                for i in range(len(data[0]) // batch_size)]
    else:
        return data


def create_vqvae(seed, input_dim, hidden_dim, z_dim, n_states):
    torch.manual_seed(seed)
    # np.random.seed(seed)
    encoder = nn.Sequential(nn.Linear(input_dim, hidden_dim),
                            nn.BatchNorm1d(hidden_dim), nn.ReLU(),
                            nn.Linear(hidden_dim, hidden_dim),
                            nn.BatchNorm1d(hidden_dim), nn.ReLU(),
                            nn.Linear(hidden_dim, z_dim))

    decoder = Decoder(input_dim, hidden_dim, z_dim)
    return VAE_VQ(encoder, decoder, z_dim, n_states)


def train_vqvae(model, training_data, n_epochs, lr, batch_size):
    optimizer = torch.optim.Adam(model.parameters(), lr)
    loss_fn = nn.GaussianNLLLoss()
    model.train()

    for i in range(n_epochs):
        for k, (feat_x, feat_y) in enumerate(batched(training_data, batch_size)):
            feat_x  = torch.tensor(feat_x, dtype=torch.float32)
            feat_y = torch.tensor(feat_y, dtype=torch.float32)

            optimizer.zero_grad()

            vq_loss, (mu, var), _, _ = model(feat_x)
            recon_loss = loss_fn(feat_y, mu, var)

            loss = recon_loss + vq_loss
            loss.backward()
            optimizer.step()

    print(f'opt_loss= {loss.item():.3}')


if __name__ == '__main__':
    # Start
    seed = int(sys.argv[1])
    data_path = sys.argv[2]
    out_dir = sys.argv[3]
    n_states = 20  # alphabet size

    # Load Data
    training_data = np.load(data_path)  # n x 10 x 2
    training_data = (training_data[:, :, 0], training_data[:, :, 1])  # (n x 10, n x 10)

    # Parameters
    input_dim = training_data[0].shape[1]
    hidden_dim = input_dim

    Z_DIM = 2
    BATCH_SIZE = 512
    LR = 1e-3
    N_EPOCHS = 4

    model =  create_vqvae(seed, input_dim, hidden_dim, Z_DIM, n_states)
    train_vqvae(model, training_data, N_EPOCHS, LR, BATCH_SIZE)
    model.eval()

    # Simplify encoder, fuse batch norms
    encoder_fused = nn.Sequential(nn.utils.fuse_conv_bn_eval(model.encoder[0], model.encoder[1]),
                                  model.encoder[2],
                                  nn.utils.fuse_conv_bn_eval(model.encoder[3], model.encoder[4]),
                                  model.encoder[5],
                                  model.encoder[6])

    # Export encoder, decoder and states
    path, name = out_dir, ''
    torch.save(encoder_fused, f'{path}/encoder{name}.pt')
    torch.save(model.decoder, f'{path}/decoder{name}.pt')
    np.savetxt(f'{path}/states{name}.txt', model.vq.embedding.weight.detach().numpy())


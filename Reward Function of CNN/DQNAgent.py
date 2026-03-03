import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import random
from collections import deque
from networks import DQN

class DQNAgent:
    def __init__(self, input_shape, num_actions, lr=1e-4, gamma=0.99, epsilon=1.0,
                 eps_min=0.1, eps_decay=1e-4, buffer_size=100000, batch_size=32):
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.num_actions = num_actions
        self.model = DQN(input_shape, num_actions).to(self.device)
        self.target_model = DQN(input_shape, num_actions).to(self.device)
        self.target_model.load_state_dict(self.model.state_dict())
        self.optimizer = optim.Adam(self.model.parameters(), lr=lr)

        self.memory = deque(maxlen=buffer_size)
        self.gamma = gamma
        self.epsilon = epsilon
        self.eps_min = eps_min
        self.eps_decay = eps_decay
        self.batch_size = batch_size
        self.step_count = 0
        self.update_freq = 1000

    def select_action(self, state):
        if np.random.rand() < self.epsilon:
            return random.randint(0, self.num_actions - 1)
        state = torch.FloatTensor(state).unsqueeze(0).to(self.device)
        with torch.no_grad():
            q_values = self.model(state)
        return torch.argmax(q_values, dim=1).item()

    def store(self, state, action, reward, next_state, done):
        self.memory.append((state, action, reward, next_state, done))
        self.epsilon = max(self.eps_min, self.epsilon - self.eps_decay)

    def train_step(self):
        if len(self.memory) < self.batch_size:
            return

        batch = random.sample(self.memory, self.batch_size)
        states, actions, rewards, next_states, dones = zip(*batch)

        states = torch.FloatTensor(np.stack(states)).to(self.device)
        next_states = torch.FloatTensor(np.stack(next_states)).to(self.device)
        actions = torch.LongTensor(actions).unsqueeze(1).to(self.device)
        rewards = torch.FloatTensor(rewards).unsqueeze(1).to(self.device)
        dones = torch.FloatTensor(dones).unsqueeze(1).to(self.device)

        q_values = self.model(states).gather(1, actions)
        with torch.no_grad():
            next_q = self.target_model(next_states).max(1)[0].unsqueeze(1)
            target = rewards + (1 - dones) * self.gamma * next_q

        loss = nn.MSELoss()(q_values, target)
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()

        self.step_count += 1
        if self.step_count % self.update_freq == 0:
            self.target_model.load_state_dict(self.model.state_dict())


import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
	b = np.load('./state_trajectory_ddp.npy')
	data_des = np.load('./state_trajectory_ddp_desired.npy')

	b = b.squeeze()
	b = b.T

	data_des = data_des.squeeze()
	data_des = data_des.T

	line_c = ['b', 'g', 'r', 'c', 'k', 'm', 'y']
	x_ = np.arange(0, b.shape[1], 1)

	# fig, axs = plt.subplots(2, 1)
	fig1 = plt.figure(1)

	for i in range(7):
		line_str = 'Joint' + str(i + 1)
		plt.plot(x_, b[i, ::], line_c[i], label=line_str)
		plt.plot(x_, data_des[i, ::], line_c[i]+'--', linewidth=2)

	plt.legend()


	fig2 = plt.figure(2)

	for i in range(3):
		line_str = 'Force' + str(i + 1)
		plt.plot(x_, b[i + 14, ::], line_c[i], label=line_str)
		plt.plot(x_, data_des[i + 14, ::], line_c[i]+'--', linewidth=2)

	plt.legend()

	# axs[0].plot(t, s1, t, s2)
	# axs[0].set_xlim(0, 2)
	# axs[0].set_xlabel('time')
	# axs[0].set_ylabel('s1 and s2')
	# axs[0].grid(True)

	# cxy, f = axs[1].cohere(s1, s2, 256, 1. / dt)
	# axs[1].set_ylabel('coherence')

	# fig.tight_layout()

	plt.show()


class hmc:

	def __init__(self, seqs, ys, sig2s, lda):
		self.seqs = seqs
		self.ys = ys
		self.sig2s = self.sig2s
		self.lda = lda
		self.A_sparse = make_A_sparse(self.seqs)
		self.E0 = dia_matrix((self.A_sparse.dot([1 / sig2s[i] for i in range(len(sig2s))]), np.array([0])), shape=(G, G))

    def hamiltonian_monte_carlo(n_samples,m,initial_position,tune=500,initial_step_size=0.1,num_steps=100,max_energy_change=1000.0,intermediate_output=True):

        num_acceptance = 0

	    initial_position = np.array(initial_position)

	    initial_potential, initial_potential_grad = potential(initial_position)

	    frac = 0.02

	    batch_size = 100

	    # collect all our samples in a list
	    samples = [initial_position]

	    # Keep a single object for momentum resampling
	    momentum = st.norm(0, m)

	    step_size = initial_step_size
	    step_size_tuning = DualAveragingStepSize(step_size)
	    # If initial_position is a 10d vector and n_samples is 100, we want 100 x 10 momentum draws
	    # we can do this in one call to np.random.normal, and iterate over rows
	    size = (n_samples + tune,) + initial_position.shape[:1]

	    for idx, p0 in enumerate(momentum.rvs(size=size)):
	        print("currently working on step " + str(idx + 1))

	        # num_steps_r = np.random.randint(int((1 - frac)*num_steps), int((1 + frac)*num_steps) + 1)
	        # step_size_r =  np.random.uniform((1 - frac)*step_size, (1 + frac)*step_size)

	        numerical_check = False

	        while numerical_check == False:

	            num_steps_r = np.random.randint(
	                int((1 - frac) * num_steps), int((1 + frac) * num_steps) + 1)
	            step_size_r = np.random.uniform(
	                (1 - frac) * step_size, (1 + frac) * step_size)

	            q_new, p_new, final_potential, final_dVdq = leapfrog(
	                samples[-1],
	                p0,
	                m,
	                potential,
	                step_size_r * num_steps_r,
	                step_size_r)

	            numerical_check = isinstance(final_potential, (float))

	        print("initial_potential = ", initial_potential)
	        print("final_potential = ", final_potential)

	        # print("log_df(p0) = ", np.sum(momentum.logpdf(p0)))
	        # print("log_df(p_new) = ", np.sum(momentum.logpdf(p_new)))

	        start_energy = -np.sum(momentum.logpdf(p0)) + initial_potential
	        new_energy = -np.sum(momentum.logpdf(p_new)) + final_potential

	        # print("new_log_p = ", str(new_log_p))

	        if new_energy == float("inf"):
	            new_energy = 1e+10

	        energy_change = new_energy - start_energy

	        if energy_change < -10:
	            energy_change = -10

	        print("energy_change = ", str(energy_change))

	        # Check Metropolis acceptance criterion
	        p_accept = min(1, np.exp(-energy_change))
	        print("p accept = ", str(np.exp(-energy_change)))
	        if np.random.rand() < p_accept:
	            samples.append(q_new)
	            initial_potential = final_potential
	            initial_potential_grad = final_dVdq
	            num_acceptance += 1
	        else:
	            samples.append(np.copy(samples[-1]))

        return np.array(samples[1 + tune:]), np.var(samples[1 + tune:], axis=0)


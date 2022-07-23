import sys


class Options:
    def __init__(self):
        self.seed = "0"
        self.N = "10"
        self.graph_density = "20"
        self.scenarios = "5"
        self.time_steps = "10"
        self.demand_range = "50:400"
        self.prod_cost_range = "200:800"
        self.unsupplied_cost = "1000"
        self.flow_cost_range = "1:1"
        self.grad_prod = "0.3"
        self.invest_flow_range = "1000:8000"
        self.invest_prod_range = "1000:8000"
        self.flow_init_max = "100"
        self.algorithm = "bilevel"

        self.unsupplied_tolerance = "1e-3"
        self.max_unsupplied = "3"

        self.cut_aggregation = "multicut"
        self.step_size = "0.5"
        self.stab_center_tol = "0.1"
        self.init_mean_value_solution = "false"

        self.heuristic_frequency = "All"
        self.heuristic_strategy = "Rand"

        self.bilevel_mode = "SOS1"
        self.primal_big_M = "1000"
        self.dual_big_M = "1000"

        self.time_limit = "21600"

    def write(self, path):
        opt_file = open(path, 'w')

        for key in vars(self):
            opt_file.write("{:<20s}= {:<20s}\n".format(key, vars(self)[key]))
        opt_file.close()


def compute_seed(N, S, T, M, max_seeds):
    for seed in range(max_seeds):
        yield(int(N) + int(S) + int(T) + int(M) + seed)


if __name__ == "__main__":

    opt_folder = sys.argv[1]

    options = Options()

    for N in ["10", "20"]:
        for S in ["1", "5", "10"]:
            for M in ["1000", "5000", "10000"]:
                for T in ["5", "10", "24"]:
                    for seed in compute_seed(N, S, T, M, 3):
                        options.N = N
                        options.seed = str(seed)
                        options.unsupplied_cost = M
                        options.scenarios = S
                        options.time_steps = T
                        if int(N) == 20:
                            options.graph_density = "10"

                        opt_file_name = "options_N" + N + "_S" + \
                            S + "_T" + T + "_M" + M + "_seed" + str(seed)
                        options.write(opt_folder + opt_file_name + ".txt")

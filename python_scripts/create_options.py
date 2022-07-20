import os


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
        self.epsilon_flow = "0.1"
        self.grad_prod = "0.1"
        self.invest_cost_range = "1000:8000"
        self.invest_prod_range = "1000:8000"
        self.flow_init_max = "100"
        self.algorithm = "bilevel"

    def write(self, path):
        opt_file = open(path, 'w')

        for key in vars(self):
            opt_file.write("{:<20s}= {:<20s}\n".format(key, vars(self)[key]))
        opt_file.close()


if __name__ == "__main__":

    opt_folder = "test_bilevel/"

    options = Options()
    for seed in range(5):
        for n in ["10", "15", "20"]:
            for scenarios in ["10", "25", "50", "100"]:
                for M in ["1000", "5000", "10000", "50000"]:
                    options.N = n
                    options.seed = str(seed)
                    options.unsupplied_cost = M
                    options.scenarios = scenarios

                    opt_file_name = "options_N" + n + "_S" + \
                        scenarios + "_seed" + str(seed) + "_M" + M
                    options.write("options/" + opt_folder +
                                  opt_file_name + ".txt")

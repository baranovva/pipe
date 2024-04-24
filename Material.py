import numpy as np
import pandas as pd


class Material:
    def __init__(self, T, p, path):
        self.T = T
        self.p = p
        self.path = "data/" + path + ".csv"
        # print(self.path)

        index_df = pd.read_csv(self.path, header=0, sep=",", usecols=[0])
        index = np.array(index_df.iloc[:, 0])

        if round(T) not in index:
            if self.path == "data/air.csv":
                material = np.array(
                    [
                        [
                            1.3037898
                            - 4.54433401e-03 * self.T
                            + 1.15712924e-05 * self.T**2
                            - 1.71098022e-08 * self.T**3
                            + 1.26728741e-11 * self.T**4
                            - 3.61429636e-15 * self.T**5,
                            1007.05722074
                            - 7.10760932e-02 * self.T
                            + 1.09243441e-03 * self.T**2
                            - 1.73622651e-06 * self.T**3
                            + 1.22716721e-09 * self.T**4
                            - 3.33705286e-13 * self.T**5,
                            0.02436655
                            + 7.87436346e-05 * self.T
                            - 1.84976615e-08 * self.T**2
                            - 2.17089760e-11 * self.T**3
                            + 1.76408367e-14 * self.T**4,
                            1.72156859e-05
                            + 5.04406359e-08 * self.T
                            - 4.16685683e-11 * self.T**2
                            + 4.90618782e-14 * self.T**3
                            - 3.82655585e-17 * self.T**4
                            + 1.21893335e-20 * self.T**5,
                            1.72049791e-05
                            + 4.94792985e-08 * self.T
                            - 3.10214433e-11 * self.T**2
                            + 1.80688009e-14 * self.T**3
                            - 4.68045260e-18 * self.T**4,
                            1/(self.T + 273)
                        ]
                    ]
                )
            elif self.path == "data/oil.csv":
                data = pd.read_csv(self.path)
                T1, T2, index1, index2 = self.find_closest_temperatures(self.T, data)
                material = np.array(
                    [
                        [
                            self.interpolate(T1, T2, Val1=data["ro"][index1], Val2=data["ro"][index2], T=T),
                            self.interpolate(T1, T2, Val1=data["c_p"][index1], Val2=data["c_p"][index2], T=T),
                            self.interpolate(T1, T2, Val1=data["lambda1"][index1], Val2=data["lambda1"][index2], T=T),
                            self.interpolate(T1, T2, Val1=data["Pr"][index1], Val2=data["Pr"][index2], T=T),
                            self.interpolate(T1, T2, Val1=data["Mu"][index1], Val2=data["Mu"][index2], T=T),
                        ]
                    ]
                )
            elif self.path == "data/methane.csv":
                data = pd.read_csv(self.path)
                T1, T2, index1, index2 = self.find_closest_temperatures(self.T, data)
                material = np.array(
                    [
                        [
                            self.interpolate(T1, T2, Val1=data["ro"][index1], Val2=data["ro"][index2], T=T),
                            self.interpolate(T1, T2, Val1=data["c_p"][index1], Val2=data["c_p"][index2], T=T),
                            self.interpolate(T1, T2, Val1=data["lambda"][index1], Val2=data["lambda"][index2], T=T),
                            self.interpolate(T1, T2, Val1=data["Pr"][index1], Val2=data["Pr"][index2], T=T),
                            self.interpolate(T1, T2, Val1=data["Mu"][index1], Val2=data["Mu"][index2], T=T),
                        ]
                    ]
                )
            else:
                data = pd.read_csv(self.path)
                T1, T2, index1, index2 = self.find_closest_temperatures(T, data)

                material = np.array(
                    [
                        [
                            999.9048951
                            + 4.79079255e-02 * self.T
                            - 7.34644522e-03 * self.T**2
                            + 3.85198135e-05 * self.T**3
                            - 1.13636364e-07 * self.T**4,
                            4210.84615385
                            - 2.03776224e00 * self.T
                            + 3.08916084e-02 * self.T**2
                            - 9.61538462e-05 * self.T**3,
                            0.54881119 
                            + 2.82137529e-03 * self.T 
                            - 1.60955711e-05 * self.T**2,
                            13.66897419
                            - 5.32406234e-01 * self.T
                            + 1.40840881e-02 * self.T**2
                            - 2.56948417e-04 * self.T**3
                            + 2.96899234e-06 * self.T**4
                            - 1.90532631e-08 * self.T**5
                            + 5.10639309e-11 * self.T**6,
                            0.00179017
                            - 6.21959139e-05 * self.T
                            + 1.80076335e-06 * self.T**2
                            - 4.56417915e-08 * self.T**3
                            + 8.57180627e-10 * self.T**4
                            - 1.00928731e-11 * self.T**5
                            + 6.44669496e-14 * self.T**6
                            - 1.68986461e-16 * self.T**7,
                            self.interpolate(T1, T2, Val1=data["beta"][index1], Val2=data["beta"][index2], T=T)]
                    ]
                )
        else:
            index = np.where(index == round(T))[0]
            material_df = pd.read_csv(self.path, header=0, sep=",", usecols=range(1, 7), skiprows=index[0], nrows=1)
            material = np.array(material_df)

        self.c_p = material[:, 1]
        self.lambd = material[:, 2]
        self.Mu = material[:, 4]
        self.Pr = material[:, 3]
        if material.shape[1] > 5:
            self.beta = material[:, 5]
        
        self.ro = material[:, 0]
        if self.path == "data/air.csv":
            self.ro *= p
            self.beta = 1/(self.T + 273)

    def interpolate(self, T1, T2, Val1, Val2, T):
        return Val1 + (Val2 - Val1) / (T2 - T1) * (T - T1)

    def find_closest_temperatures(self, T, data):
        temperatures = data["t C"].values

        closest_index = np.abs(temperatures - T).argmin()
        closest_temp = temperatures[closest_index]

        lower_index = max(closest_index - 1, 0)
        upper_index = min(closest_index + 1, len(temperatures) - 1)

        return temperatures[lower_index], temperatures[upper_index], lower_index, upper_index



if __name__ == "__main__":
    m = Material(T=0, p=1, path="water")
    print("water")
    print(f"ro {m.ro}")
    print(f"c_p {m.c_p}")
    print(f"lambda {m.lambd}")
    print(f"Pr {m.Pr}")
    print(f"Mu {m.Mu}")
    print(f"beta {m.beta}")
    m = Material(T=0, p=1, path="air")
    print("air")
    print(f"ro {m.ro}")
    print(f"c_p {m.c_p}")
    print(f"lambda {m.lambd}")
    print(f"Pr {m.Pr}")
    print(f"Mu {m.Mu}")
    print(f"beta {m.beta}")
    m = Material(T=85, p=5, path="oil")
    print("oil")
    print(f"ro {m.ro}")
    print(f"c_p {m.c_p}")
    print(f"lambda {m.lambd}")
    print(f"Pr {m.Pr}")
    print(f"Mu {m.Mu}")
    m = Material(T=85, p=5, path="methane")
    print("methane")
    print(f"ro {m.ro}")
    print(f"c_p {m.c_p}")
    print(f"lambda {m.lambd}")
    print(f"Pr {m.Pr}")
    print(f"Mu {m.Mu}")
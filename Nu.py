from functools import wraps
from buttons import show_error_popup


class Nu:
    # Среднее значение числа Нуссельта

    class NuExternal:
        # Теплоотдача при поперечном омывании цилиндрических труб. Уонг, стр 72
        def __init__(self, Re, Pr=None, is_gaz=True):
            self.Re = Re
            self.Pr = Pr
            self.is_gaz = is_gaz


        def __Nu_circle_decorator(func):
            @wraps(func)
            def wrapper(self):
                def Nu_circle(self):
                    # для цилиндров круглого сечения
                    if self.Re >= 1e-4 and self.Re <= 4e-3:
                        return 0.437 * self.Re ** 0.0895
                    elif self.Re > 4e-3 and self.Re <= 9e-2:
                        return 0.565 * self.Re ** 0.136
                    elif self.Re > 9e-2 and self.Re <= 1:
                        return 0.8 * self.Re ** 0.28
                    elif self.Re > 1 and self.Re <= 35:
                        return 0.795 * self.Re ** 0.384
                    elif self.Re > 35 and self.Re <= 5e3:
                        return 0.583 * self.Re ** 0.471
                    elif self.Re > 5e3 and self.Re <= 5e4:
                        return 0.148 * self.Re ** 0.633
                    elif self.Re > 5e4 and self.Re <= 5e5:
                        return 0.0208 * self.Re ** 0.814
                    else:
                        show_error_popup('Неверный интервал значений (Nu external)')

                return func(self) * Nu_circle(self)

            return wrapper

        @__Nu_circle_decorator
        def Nu_circle_fluid(self):
            return self.Pr ** 0.31

        @__Nu_circle_decorator
        def Nu_circle_gaz(self):
            return 1

        def calculate(self):

            if self.is_gaz:
                # если внешняя среда - газ
                return self.Nu_circle_gaz()
            else:
                # если внешняя среда - жидкость
                return 0.43 + self.Nu_circle_fluid()

        def calculate_natural(self, beta, delta_T_ln, D, nu):
            g = 9.81 
            Gr = g * beta * delta_T_ln * D**3 / nu**2
            Ra = Gr * self.Pr
            # print("beta", beta)
            # print('deltaT', delta_T_ln)
            # print("nu", nu)
            # print("Gr", Gr)
            # print("Ra", Ra)
            if Ra <= 10**9:
                return 0.8 * (Gr * self.Pr)**(1/4) * (1 + (1 + 1/(self.Pr)**0.5)**2)**(-1/4)
            else:
                phi = (self.Pr**(1/6)/(1 + 0.494 * self.Pr**(2/3)))**(2/5)
                return 0.0246 * Ra**(2/5) * phi

    class NuInternal:
        # Теплообмен при вынужденном движении в цилиндрических трубах. Уонг, стр 68
        def __init__(self, Re, Pr, is_gaz=True, Mu=None, Mu_wall=None):
            self.Re = Re
            self.Pr = Pr
            self.is_gaz = is_gaz
            self.Mu = Mu
            self.Mu_wall = Mu_wall

        def calculate(self):


            if self.is_gaz and self.Re > 2000:
                # если среда - газ
                return 0.023 * (self.Re ** 0.8) * (self.Pr ** 0.4)
            elif not self.is_gaz and (self.Pr > 0.6 and self.Pr < 100):
                # если среда - жидкость
                return 0.027 * (self.Re ** 0.8) * (self.Pr ** 0.33) * ((self.Mu / self.Mu_wall) ** 0.14)
            elif self.Re <= 2000:
                return [3.66]
            else:
                show_error_popup('Неверный интервал значений (Nu internal)')

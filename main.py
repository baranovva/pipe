import tkinter as tk
from tkinter import ttk, scrolledtext, BooleanVar, Checkbutton
from math import log, pi
from Nu import Nu
from Material import Material
from buttons import show_error_popup, show_about
from radiantion import radiation


def a(nusselt: float, lambd: float, d: float):
    return nusselt * lambd / d


def calculate():
    path_external = entry_path_external.get()
    path_internal = entry_path_internal.get()
    is_use_radiation = rad.get()
    eps = float(entry_eps.get())
    if eps < 0 or eps > 1:
        show_error_popup('Неверная степень черноты')

    if path_external == 'water':
        is_gaz_external = False
    else:
        is_gaz_external = True
    if path_internal in ['water', 'oil']:
        is_gaz_internal = False
    else:
        is_gaz_internal = True

    t_inlet = float(entry_t_inlet.get())
    t_out = float(entry_t_out.get())
    t_external = float(entry_t_external.get())

    if is_gaz_internal and (t_inlet < -50 or t_inlet > 1200):
        show_error_popup('Неверная входная температура (-50 C < T < 1200 C)')
    if is_gaz_external and (t_external < -50 or t_external > 1200):
        show_error_popup('Неверная внешняя температура (-50 C < T < 1200 C)')
    if not is_gaz_internal and (t_inlet < 0 or t_inlet > 100):
        show_error_popup('Неверная входная температура (0 C < T < 100 C)')
    if not is_gaz_external and (t_external < 0 or t_external > 100):
        show_error_popup('Неверная внешняя температура (0 C < T < 100 C)')

    if (t_out > t_inlet and t_external < t_inlet) or (t_out < t_inlet and t_external > t_inlet):
        show_error_popup('Неверная выходная температура')

    d_in = float(entry_d_in.get())
    d_external = float(entry_d_external.get())
    if d_in >= d_external or d_in == 0:
        show_error_popup('Неверные параметры трубы')

    v_external = float(entry_v_external.get())
    v_in = float(entry_v_in.get())
    if v_in <= 1e-6 or v_external <= 1e-6 or v_in > 3e+8 or v_external > 3e+8:
        show_error_popup('Неверное значение скорости')

    p_inlet = float(entry_p_inlet.get())
    p_external = float(entry_p_external.get())
    if p_inlet <= 0 or p_external <= 0:
        show_error_popup('Неверное значение давления')

    lambda_pipe = float(entry_lambda_pipe.get())
    if lambda_pipe <= 0 or lambda_pipe >= 1e+5:
        show_error_popup('Неверное значение теплопроводности трубы')

    t_avg = (t_inlet + t_out) / 2
    t_wall = (t_external + t_avg) / 2
    delta_T_max = max((t_inlet - t_external), (t_out - t_external))
    delta_T_min = min((t_inlet - t_external), (t_out - t_external))
    delta_T_ln = (delta_T_max - delta_T_min) / log(delta_T_max / delta_T_min)

    # свойства материалов при заданной температуре и давлении
    material_external = Material(T=t_external, p=p_external, path=path_external)
    material_in_wall = Material(T=t_wall, p=p_inlet, path=path_internal)
    material_in_avg = Material(T=t_avg, p=p_inlet, path=path_internal)
    material_in_inlet = Material(T=t_inlet, p=p_inlet, path=path_internal)

    Re_in = material_in_avg.ro * v_in * d_in / material_in_avg.Mu
    Re_external = material_external.ro * v_external * d_external / material_external.Mu  # Re для внешнего течения
    if use_natural_convection.get():
        if is_gaz_external:
            beta = 1/(t_external + 273)
        avg_Nu_external = Nu.NuExternal(Re=Re_external, Pr=material_external.Pr,
                                        is_gaz=is_gaz_external).calculate_natural(beta, delta_T_ln, d_external, material_external.Mu/material_external.ro)  
    else:
        avg_Nu_external = Nu.NuExternal(Re=Re_external, Pr=material_external.Pr,
                                        is_gaz=is_gaz_external).calculate()  # Уонг стр 72
    a_external_conv = a(nusselt=avg_Nu_external, lambd=material_external.lambd, d=d_external)

    a_external_rad = 0.
    if is_use_radiation and is_gaz_external:
        a_external_rad = radiation(t_wall, eps)
    a_external = a_external_conv + a_external_rad

    avg_Nu_in = Nu.NuInternal(Re=Re_in, Pr=material_in_avg.Pr, is_gaz=is_gaz_internal,
                              Mu=material_in_avg.Mu, Mu_wall=material_in_wall.Mu).calculate()  # Уонг стр 68
    a_in = a(nusselt=avg_Nu_in, lambd=material_in_avg.lambd, d=d_in)

    k_l = pi * (1 / (a_in * d_in) + log(d_external / d_in) / (2 * lambda_pipe) + 1 / (
            a_external * d_external)) ** -1  # линейный коэф теплоотдачи, Исаченко стр 37

    l = material_in_inlet.ro * v_in * material_in_inlet.c_p * pi * ((d_in / 2) ** 2) * (t_inlet - t_out) / (
            k_l * delta_T_ln)  # из уравнения теплового баланса
    delta_p = 0.184 / (Re_in ** 0.2) * (l / d_in) * (0.5 * material_in_inlet.ro * v_in ** 2)

    output_text.delete(1.0, tk.END)
    output_text.insert(tk.END, f'Re external: {Re_external[0]:.5}\n')
    output_text.insert(tk.END, f'Re internal: {Re_in[0]:.5}\n')
    output_text.insert(tk.END, f'ΔT max: {delta_T_max:.5} °C\n')
    output_text.insert(tk.END, f'ΔT min: {delta_T_min:.5} °C\n')
    output_text.insert(tk.END, f'ΔT log: {delta_T_ln:.5}\n')
    output_text.insert(tk.END, f'avarage T: {t_avg:.6} °C\n')
    output_text.insert(tk.END, f'wall T: {t_wall:.5} °C\n')
    output_text.insert(tk.END, f'Nu external: {avg_Nu_external[0]:.5}\n')
    output_text.insert(tk.END, f'Nu internal: {avg_Nu_in[0]:.5}\n')
    output_text.insert(tk.END, f'α external: {a_external[0]:.5} Вт/(м^2·K)\n')
    output_text.insert(tk.END, f'α external conv: {a_external_conv[0]:.5}  Вт/(м^2·K)\n')
    output_text.insert(tk.END, f'α external rad: {a_external_rad:.5} Вт/(м^2·K)\n')
    output_text.insert(tk.END, f'α internal: {a_in[0]:.5} Вт/(м^2·K)\n')
    output_text.insert(tk.END, f'Линейный коэффциент теплопередачи: {k_l[0]:.5} Вт/(м·K)\n')
    output_text.insert(tk.END, f'Длина трубы: {l[0]:.5} м\n')
    output_text.insert(tk.END, f'Перепад давления: {delta_p[0]:.5} Па')

    def write_to_file():
        from datetime import datetime
        file_name = str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")) + '.txt'
        with open(file_name, "w", encoding="utf-8") as file:
            file.write(f'Время расчета: {datetime.now()}\n')
            file.write(f'=====================\n')
            file.write(f'Внешняя жидкость: {path_external}\n')
            file.write(f'Внутренняя жидкость: {path_internal}\n')
            file.write(f'Внутреннее входное давление: {p_inlet} атм\n')
            file.write(f'Внешнее давление: {p_external} атм\n')
            file.write(f'=====================\n')
            file.write(f"Внутренний диаметр трубы: {d_in} м\n")
            file.write(f"Внешний диаметр трубы: {d_external} м\n")
            file.write(f"Теплопроводность трубы: {lambda_pipe} Вт/(м·K)\n")
            file.write(f'=====================\n')
            file.write(f"Скорость внешнего течения: {v_in} м/с\n")
            file.write(f"Скорость внутреннего течения: {v_external} м/с\n")
            file.write(f'Re internal: {Re_in[0]:.5}\n')
            file.write(f'Re external: {Re_external[0]:.5}\n')
            file.write(f'=====================\n')
            file.write(f'Входная температура: {t_inlet} °C\n')
            file.write(
                    f'Свойства {path_internal}: ro={material_in_inlet.ro[0]:.5} кг/м^3, c_p={material_in_inlet.c_p[0]:.5} Дж/(кг·K), lambda={material_in_inlet.lambd[0]:.5} Вт/(м·К), Pr={material_in_inlet.Pr[0]:.5}, Mu={material_in_inlet.Mu[0]:.5} H·c/м^2\n')
            file.write(f'Температура внешней среды: {t_external:.6} °C\n')
            file.write(
                    f'Свойства {path_external}: ro={material_external.ro[0]:.5} кг/м^3, c_p={material_external.c_p[0]:.5} Дж/(кг·K), lambda={material_external.lambd[0]:.5} Вт/(м·К), Pr={material_external.Pr[0]:.5}, Mu={material_external.Mu[0]:.5} H·c/м^2\n')
            file.write(f'avarage T: {t_avg:.6} °C\n')
            file.write(
                    f'Свойства {path_internal}: ro={material_in_avg.ro[0]:.5} кг/м^3, c_p={material_in_avg.c_p[0]:.5} Дж/(кг·K), lambda={material_in_avg.lambd[0]:.5} Вт/(м·К), Pr={material_in_avg.Pr[0]:.5}, Mu={material_in_avg.Mu[0]:.5} H·c/м^2\n')
            file.write(f'wall T: {t_wall:.5} °C\n')
            file.write(
                    f'Свойства {path_internal}: ro={material_in_wall.ro[0]:.5} кг/м^3, c_p={material_in_wall.c_p[0]:.5} Дж/(кг·K), lambda={material_in_wall.lambd[0]:.5} Вт/(м·К), Pr={material_in_wall.Pr[0]:.5}, Mu={material_in_wall.Mu[0]:.5} H·c/м^2\n')
            file.write(f'=====================\n')
            file.write(f'Выходная температура: {t_out} °C\n')
            file.write(f'ΔT max: {delta_T_max:.5} °C\n')
            file.write(f'ΔT min: {delta_T_min:.5} °C\n')
            file.write(f'ΔT log: {delta_T_ln:.5}\n')
            file.write(f'=====================\n')
            file.write(f'Nu external: {avg_Nu_external[0]:.5}\n')
            file.write(f'Nu internal: {avg_Nu_in[0]:.5}\n')
            file.write(f'α external: {a_external[0]:.5} Вт/(м^2·K)\n')
            file.write(f'α external conv: {a_external_conv[0]:.5}  Вт/(м^2·K)\n')
            file.write(f'α internal: {a_in[0]:.5} Вт/(м^2·K)\n')
            file.write(f'=====================\n')
            file.write(f'Учитывать радиационное изучение: {is_use_radiation}\n')
            file.write(f'Степень черноты: {eps}\n')
            file.write(f'α external rad: {a_external_rad:.5} Вт/(м^2·K)\n')
            file.write(f'=====================\n')
            file.write(f'Линейный коэффциент теплопередачи: {k_l[0]:.5} Вт/(м·K)\n')
            file.write(f'Длина трубы: {l[0]:.5} м\n')
            file.write(f'Перепад давления: {delta_p[0]:.5} Па')

    #write_to_file()


root = tk.Tk()
root.title("Конвекция 1.1.3")

label_t_inlet = tk.Label(root, text="Входная температура, °C")
label_t_inlet.grid(row=0, column=0)
entry_t_inlet = tk.Entry(root)
entry_t_inlet.grid(row=0, column=1)
entry_t_inlet.insert(0, "80")

label_t_out = tk.Label(root, text="Выходная температура, °C")
label_t_out.grid(row=1, column=0)
entry_t_out = tk.Entry(root)
entry_t_out.grid(row=1, column=1)
entry_t_out.insert(0, "50.5")

label_t_external = tk.Label(root, text="Температура внешней среды, °C")
label_t_external.grid(row=2, column=0)
entry_t_external = tk.Entry(root)
entry_t_external.grid(row=2, column=1)
entry_t_external.insert(0, "0")

label_p_inlet = tk.Label(root, text="Внутреннее входное давление, атм")
label_p_inlet.grid(row=3, column=0)
entry_p_inlet = tk.Entry(root)
entry_p_inlet.grid(row=3, column=1)
entry_p_inlet.insert(0, "1")

label_p_external = tk.Label(root, text="Внешнее давление, атм")
label_p_external.grid(row=4, column=0)
entry_p_external = tk.Entry(root)
entry_p_external.grid(row=4, column=1)
entry_p_external.insert(0, "1")

label_d_in = tk.Label(root, text="Внутренний диаметр трубы, м")
label_d_in.grid(row=5, column=0)
entry_d_in = tk.Entry(root)
entry_d_in.grid(row=5, column=1)
entry_d_in.insert(0, "0.014")

label_d_external = tk.Label(root, text="Внешний диаметр трубы, м")
label_d_external.grid(row=6, column=0)
entry_d_external = tk.Entry(root)
entry_d_external.grid(row=6, column=1)
entry_d_external.insert(0, "0.02")

label_v_external = tk.Label(root, text="Скорость внешнего течения, м/с")
label_v_external.grid(row=7, column=0)
entry_v_external = tk.Entry(root)
entry_v_external.grid(row=7, column=1)
entry_v_external.insert(0, "2")

label_v_in = tk.Label(root, text="Скорость внутреннего течения, м/с")
label_v_in.grid(row=8, column=0)
entry_v_in = tk.Entry(root)
entry_v_in.grid(row=8, column=1)
entry_v_in.insert(0, "0.089")

label_lambda_pipe = tk.Label(root, text="Теплопроводность трубы, Вт/(м·K)")
label_lambda_pipe.grid(row=9, column=0)
entry_lambda_pipe = tk.Entry(root)
entry_lambda_pipe.grid(row=9, column=1)
entry_lambda_pipe.insert(0, "0.24")

label_rad = tk.Label(root, text="Учитывать излучение")
label_rad.grid(row=10, column=0)
rad = BooleanVar()
entry_rad = Checkbutton(root, variable=rad)
entry_rad.grid(row=10, column=1)

label_eps = tk.Label(root, text="Степень черноты")
label_eps.grid(row=11, column=0)
entry_eps = tk.Entry(root)
entry_eps.grid(row=11, column=1)
entry_eps.insert(0, "0.5")

label_path_external = tk.Label(root, text="Наружная жидкость")
label_path_external.grid(row=12, column=0)
entry_path_external = ttk.Combobox(root, values=["water", "air"])
entry_path_external.grid(row=12, column=1)
entry_path_external.current(1)

label_path_internal = tk.Label(root, text="Внутренняя жидкость")
label_path_internal.grid(row=13, column=0)
entry_path_internal = ttk.Combobox(root, values=["water", "air", "oil"])
entry_path_internal.grid(row=13, column=1)
entry_path_internal.current(0)


use_natural_convection = BooleanVar()
use_natural_convection.set(False)

def toggle_external_speed():
    if use_natural_convection.get():
        entry_v_external.configure(state='disabled')
    else:
        entry_v_external.configure(state='normal')


check_use_natural_convection = Checkbutton(root, text="Использовать естественную конвекцию", variable=use_natural_convection, command=toggle_external_speed)
check_use_natural_convection.grid(row=14, columnspan=2)

button_calculate = tk.Button(root, text="Calculate", command=calculate)
button_calculate.grid(row=15, columnspan=2)

output_text = scrolledtext.ScrolledText(root, width=50, height=16, wrap=tk.WORD)
output_text.grid(row=16, columnspan=2)

button_clear = tk.Button(root, text="About", command=show_about)
button_clear.grid(row=17, columnspan=2)


root.mainloop()

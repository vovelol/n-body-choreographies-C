# n-body-choreographies-C

Этот проект реализует численное моделирование хореографий N тел в гравитационной системе. Для интегрирования системы используются два метода:
- Классический метод Рунге–Кутты 4-го порядка (RK4)
- Метод Дорманда–Принса 8-го порядка (DOPRI8)

Реализация выполнена на языке C с визуализацией 2D-анимации (используется SDL).

---

## Дифференциальное уравнение и каноническая система ОДУ

В задаче моделируется движение тел под действием гравитационного взаимодействия. Исходное дифференциальное уравнение второго порядка для каждого тела *i* имеет вид:

$$
\frac{d^2 \mathbf{r}_i}{dt^2} = \sum_{\substack{j=1 \\ j \neq i}}^{N} \frac{G M}{\|\mathbf{r}_j - \mathbf{r}_i\|^3} \, (\mathbf{r}_j - \mathbf{r}_i),
$$

где:
- \(\mathbf{r}_i = (x_i, y_i)\) — координаты тела *i*,
- \(G\) — гравитационная постоянная,
- \(M\) — масса каждого тела.

Чтобы привести систему к ОДУ первого порядка, вводим скорость

$$
\mathbf{v}_i = \frac{d \mathbf{r}_i}{dt}.
$$

Таким образом, каноническая система ОДУ, которую мы будем интегрировать, имеет следующий вид:

$$
\begin{aligned}
\frac{d \mathbf{r}_i}{dt} &= \mathbf{v}_i, \\
\frac{d \mathbf{v}_i}{dt} &= \sum_{\substack{j=1 \\ j \neq i}}^{N} \frac{G M}{\|\mathbf{r}_j - \mathbf{r}_i\|^3} \, (\mathbf{r}_j - \mathbf{r}_i),
\end{aligned}
$$

для \( i = 1, 2, \dots, N \).

---

## Физическая модель

Движение каждого тела описывается приведённой выше системой уравнений. Дополнительно:
- \( \mathbf{r}_i = (x_i, y_i) \) — координаты тела *i*.
- \( \mathbf{v}_i = (v_{x,i}, v_{y,i}) \) — скорость тела *i*.
- Уравнения учитывают гравитационное притяжение между всеми парами тел.

---

## Численные методы интегрирования

### 1. Классический метод Рунге–Кутты 4-го порядка (RK4)

Пусть \( y(t) \) — вектор состояния системы (координаты и скорости всех тел), а \( f(t, y) \) вычисляет его производную. Тогда метод RK4 задаётся формулами:

$$
\begin{aligned}
k_1 &= f(t, y), \\
k_2 &= f\!\left(t + \frac{h}{2},\; y + \frac{h}{2} k_1\right), \\
k_3 &= f\!\left(t + \frac{h}{2},\; y + \frac{h}{2} k_2\right), \\
k_4 &= f(t + h,\; y + h k_3), \\
y(t + h) &= y(t) + \frac{h}{6} \left(k_1 + 2 k_2 + 2 k_3 + k_4\right),
\end{aligned}
$$

где \( h \) — шаг интегрирования.

### 2. Метод Дорманда–Принса 8-го порядка (DOPRI8)

DOPRI8 — явный метод Рунге–Кутты с 13 этапами, обеспечивающий 8-ый порядок точности. Он вычисляет промежуточные коэффициенты следующим образом:

$$
k_i = f\!\left(t + c_i h,\; y + h \sum_{j=1}^{i-1} a_{ij} k_j\right), \quad i = 1, 2, \dots, 13,
$$

где коэффициенты \( c_i \) и \( a_{ij} \) заданы в таблицах (их конкретные значения указаны в исходном коде).

После вычисления всех \( k_i \) состояние обновляется по формуле:

$$
y(t + h) = y(t) + h \sum_{i=1}^{13} b_i k_i,
$$

где \( b_i \) — фиксированные коэффициенты.

---

## Вычисление энергии системы

Для контроля точности численного интегрирования рассчитывается полная энергия системы, состоящая из кинетической и потенциальной энергии.

### Кинетическая энергия

Для \( N \) тел кинетическая энергия вычисляется по формуле:

$$
E_{\text{kin}} = \sum_{i=1}^{N} \frac{1}{2} M \|\mathbf{v}_i\|^2
= \sum_{i=1}^{N} \frac{1}{2} M \left(v_{x,i}^2 + v_{y,i}^2\right).
$$

### Потенциальная энергия

Гравитационная потенциальная энергия системы:

$$
E_{\text{pot}} = - \sum_{i=1}^{N} \sum_{j=i+1}^{N} \frac{G M^2}{\|\mathbf{r}_j - \mathbf{r}_i\|}.
$$

### Полная энергия

Суммарная энергия системы определяется как:

$$
E = E_{\text{kin}} + E_{\text{pot}}.
$$

Контроль сохранения энергии позволяет оценить точность работы численных методов.

---

## Дополнительные детали

**Гравитационное ускорение:**

Для тела *i* ускорение, вызванное взаимодействием с телом *j*, вычисляется по формуле:

$$
\mathbf{a}_{ij} = \frac{G M}{\|\mathbf{r}_j - \mathbf{r}_i\|^3} (\mathbf{r}_j - \mathbf{r}_i).
$$

Общее ускорение для тела *i* получается суммированием вкладов от всех остальных тел.

---

**Хореографии:**

Выбранные начальные условия (например, решение "фигура-8" для трёх тел) гарантируют периодическое движение тел без столкновений, что позволяет исследовать устойчивые хореографические траектории в системе N тел.

---

## Компиляция и запуск программы

Для компиляции программы используйте следующую команду (при условии, что SDL установлена и пути указаны корректно):

```
gcc mat.c -o mat.exe -I D:\sdl\i686-w64-mingw32\include -L D:\sdl\i686-w64-mingw32\lib -lSDL3
```
Обратите внимание, что данный метод является "костыльным": пути к заголовочным файлам и библиотекам SDL указаны явно, и SDL3 необходимо установить отдельно.

После успешной компиляции запустите программу командой:
```
.\mat.exe
```

---
## Ссылки

- [Github_SDL](https://github.com/libsdl-org/SDL)
- [Library_SDL3 Download](https://github.com/libsdl-org/SDL/releases)
- [Runge–Kutta methods — Wikipedia](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
- [Dormand–Prince method — Wikipedia](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method)
- [N-body choreographies — Scholarpedia](https://web.archive.org/web/20190614005846/http://www.scholarpedia.org/article/N-body_choreographies)
- [Three-body problem — Scholarpedia](https://web.archive.org/web/20190703095653/http://www.scholarpedia.org/article/Three-body_problem)

#include <SDL3/SDL.h>
#include <stdio.h>
#include <math.h>

/*
  Программа для N=3 тел (два метода интегрирования: RK4 и DOPRI8).

  Возможные сценарии начальных условий:
    - Клавиша 0: "Лагранжев круг" (одна общая окружность для 3 тел, стандартные скорости).
    - Клавиша 1: "Лагранжев круг" с увеличенной скоростью (скорость умножается на 1.3).
    - Клавиша 2: "Наивный круг" (каждое тело на своей окружности, R=2, v = sqrt(G*M/R)).
    - Клавиша 8: "Фигура‑8".
    
  Управление отображением:
    - F1: показать оба метода (split‑screen)
    - F2: показать только RK4
    - F3: показать только DOPRI8
  
  Управление паузой:
    - F4: пауза
    - F5: продолжить
  
  Другие элементы:
    - Клавиша 7: пропустить 500 шагов (time‑skip)
    - +/-: зум, стрелки: перемещение камеры
  
  Если #define START_TIME > 0, то при выборе сценария система офлайн прокручивается на указанное время.
  При смене сценария simulation_time сбрасывается, чтобы офлайн время не накапливалось.
*/

#define START_TIME 1   // Стартовое время (сек): 0 – с начала, 6000 – с 6000 сек и т.д.
#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 600
#define N 3
#define TRAIL_LENGTH 300

typedef enum {
    MODE_BOTH,
    MODE_RK4,
    MODE_DOPRI8
} DisplayMode;
static DisplayMode currentMode = MODE_BOTH;

static const double G = 1.0;
static const double M = 1.0;
static const double h = 0.005; // шаг интегрирования

typedef struct {
    double x, y;
    double vx, vy;
} Body;

static Body bodies_RK4[N];
static Body bodies_DOPRI8[N];

static double trail_x_RK4[N][TRAIL_LENGTH];
static double trail_y_RK4[N][TRAIL_LENGTH];
static double trail_x_DOPRI8[N][TRAIL_LENGTH];
static double trail_y_DOPRI8[N][TRAIL_LENGTH];
static int trail_index = 0;

static double scale = 200.0;
static double camera_x = 0.0;
static double camera_y = 0.0;

static double simulation_time = 0.0;
static int paused = 0; // 0 – работает, 1 – пауза

// Глобальная переменная для выбора сценария (начальных условий)
static int currentScenario = 0; // Возможные значения: 0, 1, 2, 8

// Цветовые массивы для отрисовки (глобальные)
static Uint8 colorR_RK4[N] = {255, 0, 0};
static Uint8 colorG_RK4[N] = {0, 255, 0};
static Uint8 colorB_RK4[N] = {0, 0, 255};

static Uint8 colorR_DP8[N] = {255, 255, 0};
static Uint8 colorG_DP8[N] = {255, 0, 255};
static Uint8 colorB_DP8[N] = {0, 255, 255};

// Прототипы функций
static void compute_accelerations(const Body *b, double *ax, double *ay);
static void f(const Body *state, Body *dstate);
static void rk4_step(Body *state);
static void dopri8_step(Body *state);
static double compute_energy(const Body *b);
static double compute_energy_error(const Body *state1, const Body *state2);
static double compute_local_error(const Body *state1, const Body *state2);
static void reset_system(int sc);
static void world_to_screen_custom(double wx, double wy, int offsetX, int viewportWidth, int *sx, int *sy);

// ----------------------------------------------------
// Функция вычисления ускорений (без использования pow)
static void compute_accelerations(const Body *b, double *ax, double *ay) {
    for (int i = 0; i < N; i++) {
        ax[i] = 0.0;
        ay[i] = 0.0;
        for (int j = 0; j < N; j++) {
            if (j != i) {
                double dx = b[j].x - b[i].x;
                double dy = b[j].y - b[i].y;
                double r2 = dx * dx + dy * dy;
                double r = sqrt(r2);
                double r3 = r2 * r;  // Эквивалент pow(r2, 1.5)
                ax[i] += (G * M * dx) / r3;
                ay[i] += (G * M * dy) / r3;
            }
        }
    }
}

// ----------------------------------------------------
// Функция f: вычисляет производную состояния системы
// Для каждой частицы: d(state).x = state[i].vx, d(state).y = state[i].vy,
// d(state).vx = ax[i], d(state).vy = ay[i]
static void f(const Body *state, Body *dstate) {
    double ax[N], ay[N];
    compute_accelerations(state, ax, ay);
    for (int i = 0; i < N; i++) {
        dstate[i].x  = state[i].vx;
        dstate[i].y  = state[i].vy;
        dstate[i].vx = ax[i];
        dstate[i].vy = ay[i];
    }
}

// ----------------------------------------------------
// Метод Рунге–Кутты 4-го порядка (RK4)
static void rk4_step(Body *state) {
    Body k1[N], k2[N], k3[N], k4[N], temp[N];
    f(state, k1);
    for (int i = 0; i < N; i++) {
        temp[i].x  = state[i].x  + 0.5 * h * k1[i].x;
        temp[i].y  = state[i].y  + 0.5 * h * k1[i].y;
        temp[i].vx = state[i].vx + 0.5 * h * k1[i].vx;
        temp[i].vy = state[i].vy + 0.5 * h * k1[i].vy;
    }
    f(temp, k2);
    for (int i = 0; i < N; i++) {
        temp[i].x  = state[i].x  + 0.5 * h * k2[i].x;
        temp[i].y  = state[i].y  + 0.5 * h * k2[i].y;
        temp[i].vx = state[i].vx + 0.5 * h * k2[i].vx;
        temp[i].vy = state[i].vy + 0.5 * h * k2[i].vy;
    }
    f(temp, k3);
    for (int i = 0; i < N; i++) {
        temp[i].x  = state[i].x  + h * k3[i].x;
        temp[i].y  = state[i].y  + h * k3[i].y;
        temp[i].vx = state[i].vx + h * k3[i].vx;
        temp[i].vy = state[i].vy + h * k3[i].vy;
    }
    f(temp, k4);
    for (int i = 0; i < N; i++) {
        state[i].x  += (h / 6.0) * (k1[i].x + 2 * k2[i].x + 2 * k3[i].x + k4[i].x);
        state[i].y  += (h / 6.0) * (k1[i].y + 2 * k2[i].y + 2 * k3[i].y + k4[i].y);
        state[i].vx += (h / 6.0) * (k1[i].vx + 2 * k2[i].vx + 2 * k3[i].vx + k4[i].vx);
        state[i].vy += (h / 6.0) * (k1[i].vy + 2 * k2[i].vy + 2 * k3[i].vy + k4[i].vy);
    }
}

// ----------------------------------------------------
// Метод DOPRI8 (8-й порядок, фиксированный шаг)
#define S 13
static const double c[S] = {
    0.0, 1.0/18.0, 1.0/12.0, 1.0/8.0, 5.0/16.0, 3.0/8.0,
    59.0/400.0, 93.0/200.0, 5490023248.0/9719169821.0,
    13.0/20.0, 1201146811.0/1299019798.0, 1.0, 1.0
};
static const double a[S][S] = {
    {0},
    {1.0/18.0},
    {1.0/48.0, 1.0/16.0},
    {1.0/32.0, 0.0, 3.0/32.0},
    {5.0/16.0, 0.0, -75.0/64.0, 75.0/64.0},
    {3.0/80.0, 0.0, 0.0, 3.0/16.0, 3.0/20.0},
    {29443841.0/614563906.0, 0.0, 0.0, 77736538.0/692538347.0, -28693883.0/1125000000.0, 23124283.0/1800000000.0},
    {16016141.0/946692911.0, 0.0, 0.0, 61564180.0/158732637.0, 22789713.0/633445777.0, 545815736.0/2771057229.0, -180193667.0/1043307555.0},
    {39632708.0/573591083.0, 0.0, 0.0, -433636366.0/683701615.0, -421739975.0/2616292301.0, 100302831.0/723423059.0, 790204164.0/839813087.0, 800635310.0/3783071287.0},
    {246121993.0/1340847787.0, 0.0, 0.0, -37695042795.0/15268766246.0, -309121744.0/1061227803.0, -12992083.0/490766935.0, 6005943493.0/2108947869.0, 393006217.0/1396673457.0, 123872331.0/1001029789.0},
    {-1028468189.0/846180014.0, 0.0, 0.0, 8478235783.0/508512852.0, 1311729495.0/1432422823.0, -10304129995.0/1701304382.0, -48777925059.0/3047939560.0, 15336726248.0/1032824649.0, -45442868181.0/3398467696.0, 3065993473.0/597172653.0},
    {185892177.0/718116043.0, 0.0, 0.0, -3185094517.0/667107341.0, -477755414.0/1098053517.0, -703635378.0/230739211.0, 5731566787.0/1027545527.0, 5232866602.0/850066563.0, -4093664535.0/808688257.0, 3962137247.0/1805957418.0, 65686358.0/487910083.0},
    {403863854.0/491063109.0, 0.0, 0.0, -5068492393.0/434740067.0, -411421997.0/543043805.0, 652783627.0/914296604.0, 11173962825.0/925320556.0, -13158990841.0/6184727034.0, 3936647629.0/1978049680.0, -160528059.0/685178525.0, 248638103.0/1413531060.0, 0.0}
};
static const double b[S] = {
    13451932.0/455176623.0, 0.0, 0.0, 0.0, 0.0,
    -808719846.0/976000145.0, 1757004468.0/5645159321.0,
    656045339.0/265891186.0, -3867574721.0/1518517206.0,
    465885868.0/322736535.0, 53011238.0/667516719.0,
    2.0/45.0, 0.0
};

static void dopri8_step(Body *state) {
    static Body k[S][N];
    Body temp[N];
    f(state, k[0]);
    for (int i = 1; i < S; i++) {
        for (int m = 0; m < N; m++) {
            double sumx = 0.0, sumy = 0.0, sumvx = 0.0, sumvy = 0.0;
            for (int j = 0; j < i; j++) {
                sumx  += a[i][j] * k[j][m].x;
                sumy  += a[i][j] * k[j][m].y;
                sumvx += a[i][j] * k[j][m].vx;
                sumvy += a[i][j] * k[j][m].vy;
            }
            temp[m].x  = state[m].x  + h * sumx;
            temp[m].y  = state[m].y  + h * sumy;
            temp[m].vx = state[m].vx + h * sumvx;
            temp[m].vy = state[m].vy + h * sumvy;
        }
        f(temp, k[i]);
    }
    for (int m = 0; m < N; m++) {
        double sumx = 0.0, sumy = 0.0, sumvx = 0.0, sumvy = 0.0;
        for (int i = 0; i < S; i++) {
            sumx  += b[i] * k[i][m].x;
            sumy  += b[i] * k[i][m].y;
            sumvx += b[i] * k[i][m].vx;
            sumvy += b[i] * k[i][m].vy;
        }
        state[m].x  += h * sumx;
        state[m].y  += h * sumy;
        state[m].vx += h * sumvx;
        state[m].vy += h * sumvy;
    }
}

// ----------------------------------------------------
// Функция вычисления полной энергии системы (кинетическая + потенциальная)
static double compute_energy(const Body *b) {
    double kin = 0.0, pot = 0.0;
    for (int i = 0; i < N; i++) {
        double vx = b[i].vx, vy = b[i].vy;
        kin += 0.5 * (vx * vx + vy * vy);
    }
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = b[j].x - b[i].x;
            double dy = b[j].y - b[i].y;
            double r = sqrt(dx * dx + dy * dy);
            pot += -(G * M * M) / r;
        }
    }
    return kin + pot;
}

// ----------------------------------------------------
// Функция вычисления глобальной ошибки (разница энергий)
static double compute_energy_error(const Body *state1, const Body *state2) {
    double E1 = compute_energy(state1);
    double E2 = compute_energy(state2);
    return fabs(E1 - E2);
}

// ----------------------------------------------------
// Функция вычисления локальной ошибки (разница состояний)
// Вычисляется как сумма евклидовых расстояний между соответствующими координатами и скоростями для каждого тела.
static double compute_local_error(const Body *state1, const Body *state2) {
    double error = 0.0;
    for (int i = 0; i < N; i++) {
        double dx = state1[i].x - state2[i].x;
        double dy = state1[i].y - state2[i].y;
        double dvx = state1[i].vx - state2[i].vx;
        double dvy = state1[i].vy - state2[i].vy;
        error += sqrt(dx * dx + dy * dy + dvx * dvx + dvy * dvy);
    }
    return error;
}

// ----------------------------------------------------
// Функция преобразования мировых координат (wx,wy) в экранные (sx,sy)
static void world_to_screen_custom(double wx, double wy, int offsetX, int viewportWidth, int *sx, int *sy) {
    double cx = offsetX + viewportWidth * 0.5;
    double cy = WINDOW_HEIGHT * 0.5;
    double scrx = (wx - camera_x) * scale + cx;
    double scry = cy - (wy - camera_y) * scale;
    *sx = (int)scrx;
    *sy = (int)scry;
}

// ----------------------------------------------------
// Функция переустановки начальных условий (сценарий)
// Сценарий 0: Лагранжев круг (стандартный) – R = 1, v = sqrt(√3/3) ≈ 0.75984
// Сценарий 1: Лагранжев круг с увеличенной скоростью (на 30% больше)
// Сценарий 2: Наивный круг
// Сценарий 8: Фигура-8
static void reset_system(int sc) {
    currentScenario = sc;
    simulation_time = 0.0; // Сброс времени при смене сценария
    
    if (sc == 0) {
        double R = 1.0;
        double v = sqrt(sqrt(3.0) / 3.0); // ≈ 0.75984
        bodies_RK4[0].x  = R;       bodies_RK4[0].y  = 0.0;
        bodies_RK4[0].vx = 0.0;     bodies_RK4[0].vy = v;
        
        bodies_RK4[1].x  = R * cos(2.0 * M_PI / 3.0);
        bodies_RK4[1].y  = R * sin(2.0 * M_PI / 3.0);
        bodies_RK4[1].vx = -v * sin(2.0 * M_PI / 3.0);
        bodies_RK4[1].vy = v * cos(2.0 * M_PI / 3.0);
        
        bodies_RK4[2].x  = R * cos(4.0 * M_PI / 3.0);
        bodies_RK4[2].y  = R * sin(4.0 * M_PI / 3.0);
        bodies_RK4[2].vx = -v * sin(4.0 * M_PI / 3.0);
        bodies_RK4[2].vy = v * cos(4.0 * M_PI / 3.0);
    } else if (sc == 1) {
        double factor = 1.3;
        double R = 1.0;
        double v = factor * sqrt(sqrt(3.0) / 3.0);
        bodies_RK4[0].x  = R;       bodies_RK4[0].y  = 0.0;
        bodies_RK4[0].vx = 0.0;     bodies_RK4[0].vy = v;
        
        bodies_RK4[1].x  = R * cos(2.0 * M_PI / 3.0);
        bodies_RK4[1].y  = R * sin(2.0 * M_PI / 3.0);
        bodies_RK4[1].vx = -v * sin(2.0 * M_PI / 3.0);
        bodies_RK4[1].vy = v * cos(2.0 * M_PI / 3.0);
        
        bodies_RK4[2].x  = R * cos(4.0 * M_PI / 3.0);
        bodies_RK4[2].y  = R * sin(4.0 * M_PI / 3.0);
        bodies_RK4[2].vx = -v * sin(4.0 * M_PI / 3.0);
        bodies_RK4[2].vy = v * cos(4.0 * M_PI / 3.0);
    } else if (sc == 2) {
        double R = 2.0;
        double v = sqrt((G * M) / R);
        for (int i = 0; i < N; i++) {
            double theta = 2.0 * M_PI * i / N;
            bodies_RK4[i].x  = R * cos(theta);
            bodies_RK4[i].y  = R * sin(theta);
            bodies_RK4[i].vx = -v * sin(theta);
            bodies_RK4[i].vy = v * cos(theta);
        }
    } else if (sc == 8) {
        bodies_RK4[0].x  = -0.97000436;  bodies_RK4[0].y  =  0.24308753;
        bodies_RK4[0].vx =  0.466203685; bodies_RK4[0].vy =  0.43236573;
        
        bodies_RK4[1].x  =  0.97000436;  bodies_RK4[1].y  = -0.24308753;
        bodies_RK4[1].vx =  0.466203685; bodies_RK4[1].vy =  0.43236573;
        
        bodies_RK4[2].x  =  0.0;         bodies_RK4[2].y  =  0.0;
        bodies_RK4[2].vx = -0.93240737;  bodies_RK4[2].vy = -0.86473146;
    } else {
        printf("Unknown scenario %d, falling back to scenario 0.\n", sc);
        reset_system(0);
        return;
    }
    
    for (int i = 0; i < N; i++) {
        bodies_DOPRI8[i] = bodies_RK4[i];
    }
    
    if (START_TIME > 0) {
        double tskip = (double) START_TIME;
        int nsteps = (int)(tskip / h);
        printf("Offline: START_TIME = %g sec, skipping %d steps...\n", tskip, nsteps);
        for (int s = 0; s < nsteps; s++) {
            rk4_step(bodies_RK4);
            dopri8_step(bodies_DOPRI8);
            simulation_time += h;
        }
    } else {
        simulation_time = 0.0;
    }
    
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < TRAIL_LENGTH; k++) {
            trail_x_RK4[i][k] = bodies_RK4[i].x;
            trail_y_RK4[i][k] = bodies_RK4[i].y;
            trail_x_DOPRI8[i][k] = bodies_DOPRI8[i].x;
            trail_y_DOPRI8[i][k] = bodies_DOPRI8[i].y;
        }
    }
    trail_index = 0;
    
    printf("Scenario %d selected. Simulation time reset to %.3f sec.\n", currentScenario, simulation_time);
}

// ----------------------------------------------------
int main(int argc, char **argv) {
    (void)argc; (void)argv;
    
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init error: %s\n", SDL_GetError());
        return 1;
    }
    
    SDL_Window *window = SDL_CreateWindow(
        "3-Body: Scenarios (0: Lagrange, 1: Lagrange fast, 2: Naive, 8: Figure-8)",
        WINDOW_WIDTH, WINDOW_HEIGHT,
        SDL_WINDOW_RESIZABLE
    );
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
    
    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
    if (!renderer) {
        fprintf(stderr, "SDL_CreateRenderer error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    
    printf("Press one of the scenario keys: 0, 1, 2, 8\n");
    printf("Display mode: F1: Both, F2: Only RK4, F3: Only DOPRI8\n");
    printf("Pause: F4, Resume: F5\n");
    printf("Other controls: 7 => time-skip 500 steps, +/-: zoom, arrows: move camera\n");
    
    int running = 1;
    SDL_Event event;
    
    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT)
                running = 0;
            else if (event.type == SDL_EVENT_KEY_DOWN) {
                if (event.key.key == SDLK_0) {
                    reset_system(0);
                } else if (event.key.key == SDLK_1) {
                    reset_system(1);
                } else if (event.key.key == SDLK_2) {
                    reset_system(2);
                } else if (event.key.key == SDLK_8) {
                    reset_system(8);
                } else {
                    double move = 2.0 / scale;
                    switch (event.key.key) {
                        case SDLK_F1:
                            currentMode = MODE_BOTH;
                            printf("Display Mode: BOTH\n");
                            break;
                        case SDLK_F2:
                            currentMode = MODE_RK4;
                            printf("Display Mode: RK4\n");
                            break;
                        case SDLK_F3:
                            currentMode = MODE_DOPRI8;
                            printf("Display Mode: DOPRI8\n");
                            break;
                        case SDLK_F4:
                            paused = 1;
                            printf("Paused\n");
                            break;
                        case SDLK_F5:
                            paused = 0;
                            printf("Resumed\n");
                            break;
                        case SDLK_EQUALS:
                            scale *= 1.1;
                            break;
                        case SDLK_MINUS:
                            scale *= 0.9;
                            break;
                        case SDLK_UP:
                            camera_y += move;
                            break;
                        case SDLK_DOWN:
                            camera_y -= move;
                            break;
                        case SDLK_LEFT:
                            camera_x -= move;
                            break;
                        case SDLK_RIGHT:
                            camera_x += move;
                            break;
                        case SDLK_7: {
                            printf("Time-skip: 500 steps offline.\n");
                            int steps = 500;
                            for (int s = 0; s < steps; s++) {
                                rk4_step(bodies_RK4);
                                dopri8_step(bodies_DOPRI8);
                                simulation_time += h;
                            }
                            trail_index = (trail_index + 1) % TRAIL_LENGTH;
                            for (int i = 0; i < N; i++) {
                                trail_x_RK4[i][trail_index] = bodies_RK4[i].x;
                                trail_y_RK4[i][trail_index] = bodies_RK4[i].y;
                                trail_x_DOPRI8[i][trail_index] = bodies_DOPRI8[i].x;
                                trail_y_DOPRI8[i][trail_index] = bodies_DOPRI8[i].y;
                            }
                            break;
                        }
                        default:
                            break;
                    }
                }
            }
        }
        
        if (!paused) {
            if (simulation_time > 0 || currentScenario != 0) {
                rk4_step(bodies_RK4);
                dopri8_step(bodies_DOPRI8);
                simulation_time += h;
                trail_index = (trail_index + 1) % TRAIL_LENGTH;
                for (int i = 0; i < N; i++) {
                    trail_x_RK4[i][trail_index] = bodies_RK4[i].x;
                    trail_y_RK4[i][trail_index] = bodies_RK4[i].y;
                    trail_x_DOPRI8[i][trail_index] = bodies_DOPRI8[i].x;
                    trail_y_DOPRI8[i][trail_index] = bodies_DOPRI8[i].y;
                }
            }
            static double energy_timer = 0.0;
            energy_timer += h;
            if (energy_timer >= 1.0) {
                double E4 = compute_energy(bodies_RK4);
                double E8 = compute_energy(bodies_DOPRI8);
                double global_error = compute_energy_error(bodies_RK4, bodies_DOPRI8);
                double local_error = compute_local_error(bodies_RK4, bodies_DOPRI8);
                printf("t = %.3f, E_RK4 = %.6f, E_DOPRI8 = %.6f, Global_Err = %.6f, Local_Err = %.6f\n",
                       simulation_time, E4, E8, global_error, local_error);
                energy_timer = 0.0;
            }
        }
        
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        
        if (currentMode == MODE_RK4) {
            int offX = 0, vpW = WINDOW_WIDTH;
            for (int i = 0; i < N; i++) {
                SDL_SetRenderDrawColor(renderer, colorR_RK4[i], colorG_RK4[i], colorB_RK4[i], 255);
                for (int k = 0; k < TRAIL_LENGTH; k++) {
                    int idx = (trail_index + 1 + k) % TRAIL_LENGTH;
                    int sx, sy;
                    world_to_screen_custom(trail_x_RK4[i][idx], trail_y_RK4[i][idx], offX, vpW, &sx, &sy);
                    SDL_FRect fr = {(float)(sx - 1), (float)(sy - 1), 2.0f, 2.0f};
                    SDL_RenderFillRect(renderer, &fr);
                }
            }
            for (int i = 0; i < N; i++) {
                SDL_SetRenderDrawColor(renderer, colorR_RK4[i], colorG_RK4[i], colorB_RK4[i], 255);
                int sx, sy;
                world_to_screen_custom(bodies_RK4[i].x, bodies_RK4[i].y, offX, vpW, &sx, &sy);
                SDL_FRect fr = {(float)(sx - 3), (float)(sy - 3), 6.0f, 6.0f};
                SDL_RenderFillRect(renderer, &fr);
            }
        } else if (currentMode == MODE_DOPRI8) {
            int offX = 0, vpW = WINDOW_WIDTH;
            for (int i = 0; i < N; i++) {
                SDL_SetRenderDrawColor(renderer, colorR_DP8[i], colorG_DP8[i], colorB_DP8[i], 255);
                for (int k = 0; k < TRAIL_LENGTH; k++) {
                    int idx = (trail_index + 1 + k) % TRAIL_LENGTH;
                    int sx, sy;
                    world_to_screen_custom(trail_x_DOPRI8[i][idx], trail_y_DOPRI8[i][idx], offX, vpW, &sx, &sy);
                    SDL_FRect fr = {(float)(sx - 1), (float)(sy - 1), 2.0f, 2.0f};
                    SDL_RenderFillRect(renderer, &fr);
                }
            }
            for (int i = 0; i < N; i++) {
                SDL_SetRenderDrawColor(renderer, colorR_DP8[i], colorG_DP8[i], colorB_DP8[i], 255);
                int sx, sy;
                world_to_screen_custom(bodies_DOPRI8[i].x, bodies_DOPRI8[i].y, offX, vpW, &sx, &sy);
                SDL_FRect fr = {(float)(sx - 3), (float)(sy - 3), 6.0f, 6.0f};
                SDL_RenderFillRect(renderer, &fr);
            }
        } else { // MODE_BOTH (split screen)
            int halfW = WINDOW_WIDTH / 2;
            {
                int offX = 0, vpW = halfW;
                for (int i = 0; i < N; i++) {
                    SDL_SetRenderDrawColor(renderer, colorR_RK4[i], colorG_RK4[i], colorB_RK4[i], 255);
                    for (int k = 0; k < TRAIL_LENGTH; k++) {
                        int idx = (trail_index + 1 + k) % TRAIL_LENGTH;
                        int sx, sy;
                        world_to_screen_custom(trail_x_RK4[i][idx], trail_y_RK4[i][idx], offX, vpW, &sx, &sy);
                        SDL_FRect fr = {(float)(sx - 1), (float)(sy - 1), 2.0f, 2.0f};
                        SDL_RenderFillRect(renderer, &fr);
                    }
                }
                for (int i = 0; i < N; i++) {
                    SDL_SetRenderDrawColor(renderer, colorR_RK4[i], colorG_RK4[i], colorB_RK4[i], 255);
                    int sx, sy;
                    world_to_screen_custom(bodies_RK4[i].x, bodies_RK4[i].y, offX, vpW, &sx, &sy);
                    SDL_FRect fr = {(float)(sx - 3), (float)(sy - 3), 6.0f, 6.0f};
                    SDL_RenderFillRect(renderer, &fr);
                }
            }
            {
                int offX = halfW, vpW = halfW;
                for (int i = 0; i < N; i++) {
                    SDL_SetRenderDrawColor(renderer, colorR_DP8[i], colorG_DP8[i], colorB_DP8[i], 255);
                    for (int k = 0; k < TRAIL_LENGTH; k++) {
                        int idx = (trail_index + 1 + k) % TRAIL_LENGTH;
                        int sx, sy;
                        world_to_screen_custom(trail_x_DOPRI8[i][idx], trail_y_DOPRI8[i][idx], offX, vpW, &sx, &sy);
                        SDL_FRect fr = {(float)(sx - 1), (float)(sy - 1), 2.0f, 2.0f};
                        SDL_RenderFillRect(renderer, &fr);
                    }
                }
                for (int i = 0; i < N; i++) {
                    SDL_SetRenderDrawColor(renderer, colorR_DP8[i], colorG_DP8[i], colorB_DP8[i], 255);
                    int sx, sy;
                    world_to_screen_custom(bodies_DOPRI8[i].x, bodies_DOPRI8[i].y, offX, vpW, &sx, &sy);
                    SDL_FRect fr = {(float)(sx - 3), (float)(sy - 3), 6.0f, 6.0f};
                    SDL_RenderFillRect(renderer, &fr);
                }
            }
        }
        
        SDL_RenderPresent(renderer);
        SDL_Delay(16); // ~60 FPS
    }
    
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}

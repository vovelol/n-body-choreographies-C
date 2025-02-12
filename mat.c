#include <SDL3/SDL.h>
#include <stdio.h>
#include <math.h>

/*
  ОБЗОР :))))0)
  Программа для N=3 тел (двух методов: RK4 и DOPRI8).
  Есть три сценария:
    SCENARIO=0 => "круг" (одна окружность для трёх тел)
    SCENARIO=2 => "какой-то"
    SCENARIO=8 => "фигура-8"

  Дополнительно:
    #define START_TIME (сколько секунд "прокрутить" перед запуском анимации)

  Управление:
    1 => оба метода (экран делится)
    2 => только RK4
    3 => только DOPRI8
    7 => time-skip на 500 шагов
    8 => пауза
    9 => продолжить
    +/- => зум
    стрелки => сдвиг камеры
*/

#define SCENARIO 0     // 0=круг, 2=какой-то, 8=фигура-8
#define START_TIME 0   // Сколько секунд пропустить при старте (0=> начинаем с t=0, 10 => начинаем с t=10, и т.д.)

// Окно
#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 600

// N тел
#define N 3

// Длина "хвоста"
#define TRAIL_LENGTH 300 // реализация "след за точкой"

/* Режимы отображения
typedef enum {
    MODE_BOTH,
    MODE_RK4,               //РАЗОБРАТЬСЯ с SDL3 и дореализовать
    MODE_DOPRI8
} DisplayMode;

static DisplayMode currentMode = MODE_BOTH; 
*/ 

// Физика
static const double G = 1.0;
static const double M = 1.0;
static const double h = 0.005;

// Тело
typedef struct {
    double x, y;
    double vx, vy;
} Body;

static Body bodies_RK4[N];
static Body bodies_DOPRI8[N];

/*
static double trail_x_RK4[N][TRAIL_LENGTH];
static double trail_y_RK4[N][TRAIL_LENGTH];             СНАЧАЛА ДОДЕЛАТЬ ВИЗУАЛИЗАЦИЮ
static double trail_x_DOPRI8[N][TRAIL_LENGTH];
static double trail_y_DOPRI8[N][TRAIL_LENGTH];

static int trail_index = 0;
*/

/*
// Камера
static double scale    = 200.0;                                 дальность камеры(машстаб)
static double camera_x = 0.0;
static double camera_y = 0.0;                                   КОРДЫ КАМЕРЫ (0,0) всегда

// Время, пауза
static double simulation_time = 0.0;  // текущее время
static int    paused = 0;            // 0=> нет паузы, 1=> пауза
*/

// ----------------------------------------------------
// Вычисление ускорений
static void compute_accelerations(const Body *b, double *ax, double *ay) {
    for (int i=0; i<N; i++){
        ax[i] = 0.0;
        ay[i] = 0.0;
        for (int j=0;j<N;j++){
            if (j!=i){
                double dx = b[j].x - b[i].x;
                double dy = b[j].y - b[i].y;
                double r2 = dx*dx + dy*dy;
                double r3 = pow(r2,1.5);
                ax[i] += G*M*dx / r3;
                ay[i] += G*M*dy / r3;
            }
        }
    }
}
// ----------------------------------------------------
// f(state): (dx, dy, dvx, dvy)
static void f(const Body *state, Body *dstate) {
    double ax[N], ay[N];
    compute_accelerations(state, ax, ay);
    for (int i=0;i<N;i++){
        dstate[i].x  = state[i].vx;
        dstate[i].y  = state[i].vy;
        dstate[i].vx = ax[i];
        dstate[i].vy = ay[i];
    }
}
// ----------------------------------------------------
// RK4
static void rk4_step(Body *state) {
    Body k1[N], k2[N], k3[N], k4[N], temp[N];
    f(state, k1);
    for (int i=0;i<N;i++){
        temp[i].x  = state[i].x  + 0.5*h*k1[i].x;
        temp[i].y  = state[i].y  + 0.5*h*k1[i].y;
        temp[i].vx = state[i].vx + 0.5*h*k1[i].vx;
        temp[i].vy = state[i].vy + 0.5*h*k1[i].vy;
    }
    f(temp, k2);

    for (int i=0;i<N;i++){
        temp[i].x  = state[i].x  + 0.5*h*k2[i].x;
        temp[i].y  = state[i].y  + 0.5*h*k2[i].y;
        temp[i].vx = state[i].vx + 0.5*h*k2[i].vx;
        temp[i].vy = state[i].vy + 0.5*h*k2[i].vy;
    }
    f(temp, k3);

    for (int i=0;i<N;i++){
        temp[i].x  = state[i].x  + h*k3[i].x;
        temp[i].y  = state[i].y  + h*k3[i].y;
        temp[i].vx = state[i].vx + h*k3[i].vx;
        temp[i].vy = state[i].vy + h*k3[i].vy;
    }
    f(temp, k4);

    for (int i=0;i<N;i++){
        state[i].x  += (h/6.0)*(k1[i].x  + 2*k2[i].x  + 2*k3[i].x  + k4[i].x);
        state[i].y  += (h/6.0)*(k1[i].y  + 2*k2[i].y  + 2*k3[i].y  + k4[i].y);
        state[i].vx += (h/6.0)*(k1[i].vx + 2*k2[i].vx + 2*k3[i].vx + k4[i].vx);
        state[i].vy += (h/6.0)*(k1[i].vy + 2*k2[i].vy + 2*k3[i].vy + k4[i].vy);
    }
}

// ----------------------------------------------------
// DOPRI8 (8-й порядок, фиксированный шаг)
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
    for (int i=1;i<S;i++){
        for (int m=0;m<N;m++){
            double sumx=0.0, sumy=0.0, sumvx=0.0, sumvy=0.0;
            for (int j=0;j<i;j++){
                sumx  += a[i][j]*k[j][m].x;
                sumy  += a[i][j]*k[j][m].y;
                sumvx += a[i][j]*k[j][m].vx;
                sumvy += a[i][j]*k[j][m].vy;
            }
            temp[m].x  = state[m].x  + h*sumx;
            temp[m].y  = state[m].y  + h*sumy;
            temp[m].vx = state[m].vx + h*sumvx;
            temp[m].vy = state[m].vy + h*sumvy;
        }
        f(temp, k[i]);
    }

    for (int m=0;m<N;m++){
        double sumx=0.0, sumy=0.0, sumvx=0.0, sumvy=0.0;
        for (int i=0;i<S;i++){
            sumx  += b[i]*k[i][m].x;
            sumy  += b[i]*k[i][m].y;
            sumvx += b[i]*k[i][m].vx;
            sumvy += b[i]*k[i][m].vy;
        }
        state[m].x  += h*sumx;
        state[m].y  += h*sumy;
        state[m].vx += h*sumvx;
        state[m].vy += h*sumvy;
    }
}

/*
// ----------------------------------------------------
// Инициализация: scenario=0,2,8
static void init_scenario(Body *arr) {
    if (SCENARIO==0) {
        // "круг" (N=3, равносторонний треугольник)
        // double R=2.0;
        // double v=sqrt(G*M/(3*R));
        // arr[0].x=..., arr[0].vx=..., ...
        // ...
        printf("SCENARIO=0 => circle (one circle for 3 bodies)\n");
    }
    else if (SCENARIO==2) {
        какой-то
        printf("SCENARIO=2 => chto-to da budet \n");
    }
    else if (SCENARIO==8) {
        // "Фигура-8" (three-body figure-8)
        // arr[0].x=..., arr[0].vx=..., ...
    }
}
*/

/*
// ----------------------------------------------------
// Преобразование (wx,wy)->(sx,sy)
static void world_to_screen_custom(
    ЧЕКНУТЬ ВИДОС с разбором
){
   типо реализация из кордов в экран
}
*/

// ----------------------------------------------------
// Энергия
static double compute_energy(const Body *b) {
    double kin=0.0, pot=0.0;
    for (int i=0;i<N;i++){
        double vx=b[i].vx, vy=b[i].vy;
        kin += 0.5*(vx*vx+vy*vy);
    }
    for (int i=0;i<N;i++){
        for (int j=i+1;j<N;j++){
            double dx=b[j].x-b[i].x;
            double dy=b[j].y-b[i].y;
            double r=sqrt(dx*dx+dy*dy);
            pot += -(G*M*M)/r;
        }
    }
    return kin+pot;
}

// ----------------------------------------------------
// Цвета
static Uint8 colorR_RK4[N] = {255,   0,   0};
static Uint8 colorG_RK4[N] = {  0, 255,   0};
static Uint8 colorB_RK4[N] = {  0,   0, 255};

static Uint8 colorR_DP8[N] = {255, 255,   0};
static Uint8 colorG_DP8[N] = {255,   0, 255};
static Uint8 colorB_DP8[N] = {  0, 255, 255};

// ----------------------------------------------------
int main(int argc, char **argv){
    (void)argc; (void)argv;

   /*
    типо вызываем SDL
    реализуем
    окошко
    туда сюда
   
   */

    // 1) Инициализируем начальное состояние
    init_scenario(bodies_RK4);
    init_scenario(bodies_DOPRI8);

    /*
    // 1.1) "START_TIME" офлайновая прокрутка
#if (START_TIME > 0)
    {
        типо скипаем по кадрам в зависимости от h
        реализовать в функциях её
    }
#endif
    */

    /*
    // 2) Заполняем хвост (начальное положение)
    for (int i=0;i<N;i++){
        туда сюда делаем хвост-след 
        для двух систем
        i k = боди x
        i k = боди y
        т.д.
    }
    trail_index=0;
    */

    double energy_timer=0.0;
    int running=1;
    SDL_Event event;

    printf("SCENARIO=%d, START_TIME=%d\n", SCENARIO, START_TIME);
    printf("Controls:\n");
    printf("1 => both, 2 => only RK4, 3 => only DOPRI8\n");
    printf("7 => time-skip 500 steps\n");
    printf("8 => pause, 9 => resume\n");
    printf("+/- => zoom, arrows => move camera\n");

    /*
    while (running){
        // События +- 200 строк :(((

        Много кейсой 0_о ужс
        зум
        вызов окна и смена на 1(2),2(RK),3(dopri)

        в консоль импульс писать за некоторые сек.

        в текстовый файл 2 отдельных писать импульс

        7(skip)
        8(zoom)
        9(zoom)
        ужасный рендер и их всего

        нужен хвост, положение, синхра
    }
    */
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}

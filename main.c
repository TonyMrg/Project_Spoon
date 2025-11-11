#include <stdio.h>
#include <math.h>
#include <SDL.h>
#include <stdbool.h>
#define N_STEPS 20000

void find_accelerations(double a, double b, double c, double d, double c1, double c2, double* th1, double* th2);
void play_double_pendulum(const double* t,
    const double* theta1,
    const double* theta2,
    int n,
    double L1,
    double L2,
    int width,
    int height);
static void draw_double_pendulum(SDL_Renderer* renderer,
    double theta1, double theta2,
    double L1, double L2,
    int width, int height)
{
    double scale = (height / 2.5) / (L1 + L2);

    int x0 = width / 2;
    int y0 = height / 4;

    int x1 = x0 + (int)(L1 * scale * sin(theta1));
    int y1 = y0 + (int)(L1 * scale * cos(theta1));

    int x2 = x1 + (int)(L2 * scale * sin(theta2));
    int y2 = y1 + (int)(L2 * scale * cos(theta2));

    // Clear screen
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);

    // Draw rods
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderDrawLine(renderer, x0, y0, x1, y1);
    SDL_RenderDrawLine(renderer, x1, y1, x2, y2);

    // Draw bobs (simple filled circles)
    int radius = 8;
    SDL_SetRenderDrawColor(renderer, 200, 0, 0, 255);
    for (int w = -radius; w <= radius; w++) {
        for (int h = -radius; h <= radius; h++) {
            if (w * w + h * h <= radius * radius)
                SDL_RenderDrawPoint(renderer, x1 + w, y1 + h);
        }
    }
    SDL_SetRenderDrawColor(renderer, 0, 0, 200, 255);
    for (int w = -radius; w <= radius; w++) {
        for (int h = -radius; h <= radius; h++) {
            if (w * w + h * h <= radius * radius)
                SDL_RenderDrawPoint(renderer, x2 + w, y2 + h);
        }
    }

    SDL_RenderPresent(renderer);
}

int main() {
    double m1 = 1;
    double m2 = 1;
    double l1 = 1;
    double l2 = 1;
    double theta1 = 1.5708;
    double theta2 = 1.5708;
    double theta1t = 0.1;
    double theta2t = 0.1;

    double theta1_arr[N_STEPS];
    double theta2_arr[N_STEPS];
    double t_arr[N_STEPS];

    double t = 0.0;
    for (int i = 0; i < N_STEPS; i++, t += 0.01) {
        double a = (m1 + m2) * l1;
        double b = m2 * l2 * cos(theta1 - theta2);
        double c = m2 * l1 * cos(theta1 - theta2);
        double d = m2 * l2;
        double c1 = m2 * l2 * pow(theta2t, 2) * sin(theta1 - theta2) + (m1 + m2) * 9.81 * sin(theta1);
        double c2 = -m2 * l1 * pow(theta1t, 2) * sin(theta1 - theta2) + m2 * 9.81 * sin(theta2);
        double theta1tt;
        double theta2tt;
        find_accelerations(a, b, c, d, c1, c2, &theta1tt, &theta2tt);

        theta1 = theta1 + theta1t * 0.01;
        theta2 = theta2 + theta2t * 0.01;
        theta1t = theta1t + theta1tt * 0.01;
        theta2t = theta2t + theta2tt * 0.01;

        theta1_arr[i] = theta1;
        theta2_arr[i] = theta2;
        t_arr[i] = t;
    }

    // Example usage: play animation
    play_double_pendulum(t_arr, theta1_arr, theta2_arr, N_STEPS, l1, l2, 800, 600);

    return 0;
}
void find_accelerations(double a, double b, double c, double d, double c1, double c2, double* th1, double* th2) {
    *th1 = 1 / (a * d - b * c) * (-d * c1 + b * c2);
    *th2 = 1 / (a * d - b * c) * (c * c1 - a * c2);
}

void play_double_pendulum(const double* t,
    const double* theta1,
    const double* theta2,
    int n,
    double L1,
    double L2,
    int width,
    int height)
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("Double Pendulum",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        width, height, 0);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1,
        SDL_RENDERER_ACCELERATED);

    bool running = true;
    SDL_Event e;
    Uint32 start_time = SDL_GetTicks();

    for (int i = 0; i < n && running; i++) {
        // Handle events
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = false;
        }

        draw_double_pendulum(renderer, theta1[i], theta2[i], L1, L2, width, height);

        if (i < n - 1) {
            double dt = t[i + 1] - t[i];  // seconds
            SDL_Delay((Uint32)(dt * 1000)); // convert to ms
        }
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}
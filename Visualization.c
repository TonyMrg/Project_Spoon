#include "global_parameters.h"


void visualize_double_pendulum(const double* results, int steps, double dt, double playback_time_s, int n_masses)
{
    // --- Window settings (adjust if desired) ---
    const int WIDTH = 2400;    // set your window width
    const int HEIGHT = 1300;   // set your window height

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL init failed: %s\n", SDL_GetError());
        return;
    }

    SDL_Window* window = SDL_CreateWindow(
        "Pendulum Visualization",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIDTH, HEIGHT, SDL_WINDOW_SHOWN);

    if (!window) {
        fprintf(stderr, "Window creation failed: %s\n", SDL_GetError());
        SDL_Quit();
        return;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer) {
        fprintf(stderr, "Renderer creation failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return;
    }

    // ----------------------------------------------------------
    // Dynamic scaling to keep pendulum on screen
    // ----------------------------------------------------------
    double total_length_m = n_masses * 1.0; // if each link is 1 meter long
    double max_screen_height_px = HEIGHT * 0.8;
    double SCALE = max_screen_height_px / total_length_m;

    const int CX = WIDTH / 2;
    const int CY = HEIGHT / 6; // lower the pivot a bit for visibility
    // ----------------------------------------------------------

    int running = 1;
    SDL_Event e;

    // ----------------------------------------------------------
    // Playback settings (adjust playback_time_s or desired_fps)
    // ----------------------------------------------------------
    double simulated_duration_s = steps * dt;
    double desired_fps = 165.0; // adjust playback frame rate if needed
    int total_frames = (int)(playback_time_s * desired_fps);

    if (total_frames < 1) total_frames = 1;

    double step_increment = (double)steps / (double)total_frames;
    if (step_increment < 1.0) step_increment = 1.0;

    int delay_ms = (int)(1000.0 / desired_fps);
    // ----------------------------------------------------------

    SDL_Color mass_colors[] = {
        {255, 80, 80, 255}, {80, 160, 255, 255}, {80, 255, 120, 255},
        {255, 200, 50, 255}, {180, 100, 255, 255}, {255, 100, 180, 255}
    };
    int n_colors = sizeof(mass_colors) / sizeof(SDL_Color);

    for (double f = 0; f < steps && running; f += step_increment) {
        int step = (int)f;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = 0;
        }

        SDL_SetRenderDrawColor(renderer, 20, 20, 30, 255);
        SDL_RenderClear(renderer);

        int prev_px = CX, prev_py = CY;
        for (int i = 0; i < n_masses; ++i) {
            double x = results[step * n_masses * 2 + 2 * i];
            double y = results[step * n_masses * 2 + 2 * i + 1];

            int px = CX + (int)(SCALE * x);
            int py = CY - (int)(SCALE * y);

            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL_RenderDrawLine(renderer, prev_px, prev_py, px, py);

            SDL_Rect mrect = { px, py, 3, 3 };
            SDL_Color c = mass_colors[i % n_colors];
            SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, c.a);
            SDL_RenderFillRect(renderer, &mrect);

            prev_px = px;
            prev_py = py;
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(delay_ms);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

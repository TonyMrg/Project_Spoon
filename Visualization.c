#include "global_parameters.h"
#define MAX_TRAIL 200 // number of past positions to store for each mass

void visualize_double_pendulum(const double* results, int steps, double dt, double playback_time_s, int n_masses)
{
    // --- Window / GL settings ---
    const int WIDTH = 1600;
    const int HEIGHT = 1000;

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL init failed: %s\n", SDL_GetError());
        return;
    }

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

    SDL_Window* window = SDL_CreateWindow(
        "Pendulum 3D Visualization",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIDTH, HEIGHT,
        SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);

    if (!window) {
        fprintf(stderr, "Window creation failed: %s\n", SDL_GetError());
        SDL_Quit();
        return;
    }

    SDL_GLContext glctx = SDL_GL_CreateContext(window);
    if (!glctx) {
        fprintf(stderr, "Failed to create GL context: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return;
    }

    glViewport(0, 0, WIDTH, HEIGHT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_NORMALIZE);
    glClearColor(0.078f, 0.078f, 0.118f, 1.0f);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    {
        GLfloat light_pos[] = { 1.0f, 1.0f, 2.0f, 0.0f };
        GLfloat ambient[] = { 0.10f, 0.10f, 0.10f, 1.0f };
        GLfloat diffuse[] = { 0.9f, 0.9f, 0.9f, 1.0f };
        glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    }

    SDL_Color mass_colors[] = {
        {255, 80, 80, 255}, {80, 160, 255, 255}, {80, 255, 120, 255},
        {255, 200, 50, 255}, {180, 100, 255, 255}, {255, 100, 180, 255}
    };
    int n_colors = sizeof(mass_colors) / sizeof(SDL_Color);

    double total_length_m = n_masses * 1.0;
    double SCALE = 1.0;

    double simulated_duration_s = steps * dt;
    double desired_fps = 165.0;
    int total_frames = (int)(playback_time_s * desired_fps);
    if (total_frames < 1) total_frames = 1;
    double step_increment = (double)steps / (double)total_frames;
    if (step_increment < 1.0) step_increment = 1.0;
    int delay_ms = (int)(1000.0 / desired_fps);

    int running = 1;
    SDL_Event e;

    GLUquadric* quad = gluNewQuadric();

    // Variables for rotation and zoom
    double cam_dist = total_length_m * 3.0 * SCALE; // Initial zoom level
    double cam_angle_x = 0.0; // Horizontal rotation angle
    double cam_angle_y = 0.0; // Vertical rotation angle
    int rotating = 0;         // Is the left mouse button pressed?
    int last_mouse_x = 0, last_mouse_y = 0;

    while (running) {
        for (double f = 0; f < steps && running; f += step_increment) {
            int step = (int)f;

            while (SDL_PollEvent(&e)) {
                if (e.type == SDL_QUIT) {
                    running = 0;
                }
                else if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_RESIZED) {
                    int w = e.window.data1;
                    int h = e.window.data2;
                    if (h == 0) h = 1;
                    glViewport(0, 0, w, h);
                }
                else if (e.type == SDL_KEYDOWN) {
                    if (e.key.keysym.sym == SDLK_ESCAPE) running = 0;
                }
                else if (e.type == SDL_MOUSEBUTTONDOWN) {
                    if (e.button.button == SDL_BUTTON_LEFT) {
                        rotating = 1;
                        last_mouse_x = e.button.x;
                        last_mouse_y = e.button.y;
                    }
                }
                else if (e.type == SDL_MOUSEBUTTONUP) {
                    if (e.button.button == SDL_BUTTON_LEFT) {
                        rotating = 0;
                    }
                }
                else if (e.type == SDL_MOUSEMOTION) {
                    if (rotating) {
                        int dx = e.motion.x - last_mouse_x;
                        int dy = e.motion.y - last_mouse_y;
                        cam_angle_x -= dx * 0.2; // Inverted horizontal rotation
                        cam_angle_y += dy * 0.2;
                        if (cam_angle_y > 89.0) cam_angle_y = 89.0; // Limit vertical rotation
                        if (cam_angle_y < -89.0) cam_angle_y = -89.0;
                        last_mouse_x = e.motion.x;
                        last_mouse_y = e.motion.y;
                    }
                }
                else if (e.type == SDL_MOUSEWHEEL) {
                    cam_dist -= e.wheel.y * 5; // Adjust zoom sensitivity
                    if (cam_dist < 1.0) cam_dist = 1.0; // Prevent zooming too close
                }
            }

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            int w, h;
            SDL_GetWindowSize(window, &w, &h);
            float aspect = (h == 0) ? 1.0f : (float)w / (float)h;
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluPerspective(45.0, aspect, 0.01, 1000.0);

            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();

            // Calculate camera position based on rotation and zoom
            double cam_x = cam_dist * cos(cam_angle_y * M_PI / 180.0) * sin(cam_angle_x * M_PI / 180.0);
            double cam_z = cam_dist * cos(cam_angle_y * M_PI / 180.0) * cos(cam_angle_x * M_PI / 180.0);
            double cam_y = cam_dist * sin(cam_angle_y * M_PI / 180.0);

            gluLookAt(cam_x, cam_y, cam_z,
                      0.0, 0.0, 0.0,
                      0.0, 1.0, 0.0);

            // Draw the floor grid
            glDisable(GL_LIGHTING);
            glColor3f(0.5f, 0.5f, 0.5f); // Gray color for the grid
            glLineWidth(1.0f);
            glBegin(GL_LINES);
            float grid_y = (float)(-n_masses); // Grid at y = -n_masses
            int grid_size = 10; // Number of lines in each direction
            float grid_spacing = (float)n_masses; // Scale spacing with n_masses
            for (int i = -grid_size; i <= grid_size; ++i) {
                float pos = i * grid_spacing;
                glVertex3f(pos, grid_y, -grid_size * grid_spacing);
                glVertex3f(pos, grid_y, grid_size * grid_spacing);
                glVertex3f(-grid_size * grid_spacing, grid_y, pos);
                glVertex3f(grid_size * grid_spacing, grid_y, pos);
            }
            glEnd();
            glEnable(GL_LIGHTING);

            // Draw the pendulum (same as before)
            float prev_x = 0.0f, prev_y = 0.0f, prev_z = 0.0f;
            glLineWidth(2.0f);
            glDisable(GL_LIGHTING);
            glBegin(GL_LINES);
            for (int i = 0; i < n_masses; ++i) {
                double x = results[step * n_masses * 3 + 3 * i + 0];
                double y = results[step * n_masses * 3 + 3 * i + 1];
                double z = results[step * n_masses * 3 + 3 * i + 2];

                float fx = (float)(SCALE * x);
                float fy = (float)(SCALE * y);
                float fz = (float)(SCALE * z);

                glVertex3f(prev_x, prev_y, prev_z);
                glVertex3f(fx, fy, fz);

                prev_x = fx; prev_y = fy; prev_z = fz;
            }
            glEnd();
            glEnable(GL_LIGHTING);

            prev_x = 0.0f; prev_y = 0.0f; prev_z = 0.0f;
            for (int i = 0; i < n_masses; ++i) {
                double x = results[step * n_masses * 3 + 3 * i + 0];
                double y = results[step * n_masses * 3 + 3 * i + 1];
                double z = results[step * n_masses * 3 + 3 * i + 2];

                float fx = (float)(SCALE * x);
                float fy = (float)(SCALE * y);
                float fz = (float)(SCALE * z);

                SDL_Color c = mass_colors[i % n_colors];
                GLfloat mat_diffuse[4] = { c.r / 255.0f, c.g / 255.0f, c.b / 255.0f, 1.0f };
                GLfloat mat_specular[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
                GLfloat mat_shininess[1] = { 20.0f };
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

                glPushMatrix();
                glTranslatef(fx, fy, fz);
                double radius = 0.06 * SCALE;
                gluSphere(quad, (GLdouble)radius, 16, 12);
                glPopMatrix();

                prev_x = fx; prev_y = fy; prev_z = fz;
            }

            glDisable(GL_LIGHTING);
            glPointSize(6.0f);
            glBegin(GL_POINTS);
            glColor3f(1.0f, 1.0f, 1.0f);
            glVertex3f(0.0f, 0.0f, 0.0f);
            glEnd();
            glEnable(GL_LIGHTING);

            SDL_GL_SwapWindow(window);
            SDL_Delay(delay_ms);
        }
        running = 0;
    }

    gluDeleteQuadric(quad);
    SDL_GL_DeleteContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
}
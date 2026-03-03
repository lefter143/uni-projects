/* ELEFTHERIOS PANTZARTZIS 4344 */
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>



typedef GLfloat point2[2]; /* define a point data type */
point2 points[20000];
GLfloat pointColors[20000][3];
int numPoints = 15000, lastMouseX, lastMouseY, multiplier = 1, isDragging = 0, randomColor = 0;
float newX = 0.0f, newY = 0.0f;

void myInit(void) {
    /* attributes */
    glEnable(GL_BLEND);
    glClearColor(1.0, 1.0, 1.0, 0.0); /* white background */
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT);  /* clear the window */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-5.0 * multiplier + newX, 5.0 * multiplier + newX,
        multiplier + newY, 10.0 * multiplier + newY); /* smaller bounds because the fern was too small and wrongly sized */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glBegin(GL_POINTS);
    for (int t = 0; t < numPoints; t++) {
        glColor3f(pointColors[t][0], pointColors[t][1], pointColors[t][2]);
        glVertex2fv(points[t]);
    }
    glEnd();
    glutSwapBuffers();
}

void drawFernPoints() {
    float x = 0.0f, y = 0.0f, xn = 0.0f, yn = 0.0f;
    /* computes and plots numPoints new points */
    for (int t = 0; t < numPoints; t++) {
        float r = (float)rand() / RAND_MAX;

        if (r < 0.01) {
            xn = 0.0;
            yn = 0.16 * y;
        }
        else if (r < 0.86) {
            xn = 0.85 * x + 0.04 * y;
            yn = -0.04 * x + 0.85 * y + 1.6;
        }
        else if (r < 0.93) {
            xn = 0.2 * x - 0.26 * y;
            yn = 0.23 * x + 0.22 * y + 1.6;
        }
        else {
            xn = -0.15 * x + 0.28 * y;
            yn = 0.26 * x + 0.24 * y + 0.44;
        }
        x = xn;
        y = yn;
        points[t][0] = x;
        points[t][1] = y;

        if (randomColor) {
            pointColors[t][0] = (float)(rand() % 256) / 255.0f;
            pointColors[t][1] = (float)(rand() % 256) / 255.0f;
            pointColors[t][2] = (float)(rand() % 256) / 255.0f;
        }
        else {
            pointColors[t][0] = 1.0f;
            pointColors[t][1] = 0.0f;
            pointColors[t][2] = 0.0f;
        }
    }
}

void menu(int option) {
    switch (option) {
    case 1:
        numPoints = 20000;
        randomColor = 1;
        break;
    case 2:
        numPoints = 15000;
        randomColor = 0;
        break;
    case 3:
        multiplier *= 2;
        break;
    case 4:
        multiplier = 1;
        break;
    case 5:
        exit(0);
    }
    drawFernPoints();
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            isDragging = 1;
            lastMouseX = x;
            lastMouseY = y;
        }
        else {
            isDragging = 0;
        }
    }
}

void motion(int x, int y) {
    if (isDragging) {
        float dx = (x - lastMouseX) * 0.02f;
        float dy = (y - lastMouseY) * 0.02f;

        newX -= dx;
        newY += dy;

        lastMouseX = x;
        lastMouseY = y;
        glutPostRedisplay();
    }
}

int main(int argc, char** argv) {
    /* Standard GLUT initialization */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(500, 500); /* 500 x 500 pixel window */
    glutInitWindowPosition(0, 0); /* place window top left on display */
    glutCreateWindow("Lorem Ipsum"); /* window title */
    myInit(); /* set attributes */

    drawFernPoints();
    glutDisplayFunc(display); /* display callback invoked when window opened */
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutCreateMenu(menu);
    glutAddMenuEntry("20000 points, random color", 1);
    glutAddMenuEntry("15000 points, same color", 2);
    glutAddMenuEntry("Cut size in half", 3);
    glutAddMenuEntry("Initial size", 4);
    glutAddMenuEntry("Finish", 5);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glutMainLoop(); /* enter event loop */
}

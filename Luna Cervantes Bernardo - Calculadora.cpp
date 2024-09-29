#define _USE_MATH_DEFINES

// #define M_PI 3.14159265358979323846264338327950288419716939937510

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <complex>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <functional>
#include <limits>

typedef std::function<double(double)> FuncionReal;
typedef std::function<double(std::function<double(double)>, double)> DerivadaFR;

//---------------------------- metodos de uso general ------------------------------------

void imprimirPolinomio(int grado, const std::vector<double> &coeficientes)
{

    std::cout << std::fixed << std::setprecision(2);

    for (int i = grado; i >= 0; i--)
    {
        if (coeficientes[grado - i] != 0.0)
        {
            if (i == grado)
            {
                if (coeficientes[grado - i] < 0)
                {
                    std::cout << "-";
                }
                std::cout << std::abs(coeficientes[grado - i]) << "x^" << i;
            }
            else
            {
                if (coeficientes[grado - i] < 0)
                {
                    std::cout << " - ";
                }
                else
                {
                    std::cout << " + ";
                }
                std::cout << std::abs(coeficientes[grado - i]) << "x^" << i;
            }
        }
    }
}

int obtenerSigno(double numero)
{
    if (numero > 0)
    {
        return 1;
    }
    else if (numero < 0)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

double evaluarPolinomio(const std::vector<double> &coeficientes, double x)
{
    double resultado = 0.0;
    for (size_t i = 0; i < coeficientes.size(); ++i)
    {
        resultado += coeficientes[i] * std::pow(x, coeficientes.size() - 1 - i);
    }
    return resultado;
}

double evaluarDerivada(const std::vector<double> &coeficientes, double x)
{
    double resultado = 0.0;
    for (size_t i = 0; i < coeficientes.size() - 1; ++i)
    {
        // n*ax^(n-1)
        resultado += (coeficientes.size() - 1 - i) * coeficientes[i] * std::pow(x, coeficientes.size() - 2 - i);
    }
    return resultado;
}

bool esCuadradoPerfecto(double a, double b, double c, double tolerancia = 1e-6)
{
    // Un polinomio cuadrático ax^2 + bx + c es un cuadrado perfecto si b^2 = 4ac
    return std::abs(b * b - 4 * a * c) < tolerancia;
}

double atan2EnGrados(double y, double x)
{
    double radianes = atan2(y, x);
    return radianes * 180.0 / M_PI;
}

double senoEnGrados(double grados)
{
    double radianes = grados * M_PI / 180.0;
    return sin(radianes);
}

double cosenoEnGrados(double grados)
{
    double radianes = grados * M_PI / 180.0;
    return cos(radianes);
}

std::pair<std::complex<double>, std::complex<double>> resolverx_2(double a, double b, double c)
{
    std::complex<double> raiz1, raiz2;
    double discriminante = (b * b) - (4 * a * c);

    if (discriminante >= 0)
    {
        raiz1 = std::complex<double>((-b + sqrt(discriminante)) / (2 * a), 0);
        raiz2 = std::complex<double>((-b - sqrt(discriminante)) / (2 * a), 0);
    }
    else
    {
        double parteReal = -b / (2 * a);
        double parteImaginaria = sqrt(-discriminante) / (2 * a);
        raiz1 = std::complex<double>(parteReal, parteImaginaria);
        raiz2 = std::complex<double>(parteReal, -parteImaginaria);
    }

    return std::make_pair(raiz1, raiz2);
}

double resolverx_1(double a, double b)
{
    double raiz = (-b) / a;
    return raiz;
}

//------------------------------ division sintetica --------------------------------------

std::vector<double> divisionSintetica(double a, double b, double c, double d, double raiz) // grado 3
{

    std::vector<double> nuevoPolinomio(3);

    /*

    std::cout << "a: " << a << "\n";
    std::cout << "b: " << b << "\n";
    std::cout << "c: " << c << "\n";
    std::cout << "d: " << d << "\n";
    std::cout << "raiz: " << raiz << "\n";

    */

    nuevoPolinomio[0] = a;
    nuevoPolinomio[1] = (b) + (nuevoPolinomio[0] * raiz);
    nuevoPolinomio[2] = (c) + (nuevoPolinomio[1] * raiz);

    /*

    std::cout << "nuevoPolinomio[0]: " << nuevoPolinomio[0] << "\n";
    std::cout << "nuevoPolinomio[1]: " << nuevoPolinomio[1] << "\n";
    std::cout << "nuevoPolinomio[2]: " << nuevoPolinomio[2] << "\n";

    */

    double comprobacion = d + (nuevoPolinomio[2] * raiz);

    // std::cout << "La operacion final que se realizo dio como resultado: " << comprobacion << "\n";

    if (std::abs(comprobacion) <= 1e-10)
    {
        return nuevoPolinomio;
    }
    else
    {
        std::cout << "La division sintetica no se aproxima a cero, aunque los calculos podrian funcionar\n";
        std::cout << "El residuo es : " << comprobacion << "\n";
        std::cout << "\n";
        return nuevoPolinomio;
    }
}

std::vector<double> divisionSintetica(double a, double b, double c, double d, double e, double raiz) // grado 4
{

    std::vector<double> nuevoPolinomio(4);

    /*

    std::cout << "a: " << a << "\n";
    std::cout << "b: " << b << "\n";
    std::cout << "c: " << c << "\n";
    std::cout << "d: " << d << "\n";
    std::cout << "raiz: " << raiz << "\n";

    */

    nuevoPolinomio[0] = a;
    nuevoPolinomio[1] = (b) + (nuevoPolinomio[0] * raiz);
    nuevoPolinomio[2] = (c) + (nuevoPolinomio[1] * raiz);
    nuevoPolinomio[3] = (d) + (nuevoPolinomio[2] * raiz);

    /*

    std::cout << "nuevoPolinomio[0]: " << nuevoPolinomio[0] << "\n";
    std::cout << "nuevoPolinomio[1]: " << nuevoPolinomio[1] << "\n";
    std::cout << "nuevoPolinomio[2]: " << nuevoPolinomio[2] << "\n";
    std::cout << "nuevoPolinomio[2]: " << nuevoPolinomio[3] << "\n";

    */

    double comprobacion = e + (nuevoPolinomio[3] * raiz);

    // std::cout << "La operacion final que se realizo dio como resultado: " << comprobacion << "\n";

    if (std::abs(comprobacion) <= 1e-10)
    {
        return nuevoPolinomio;
    }
    else
    {
        std::cout << "La division sintetica no se aproxima a cero, aunque los calculos podrian funcionar\n";
        std::cout << "El residuo es : " << comprobacion << "\n";
        std::cout << "\n";
        return nuevoPolinomio;
    }
}

//---------------------------- division de polinomios ------------------------------------

// las divisiones solo funcionan si se esta dividiendo entre un polinomio de grado 2

std::vector<double> dividirPolinomios(double a, double b, double r) // division grado 3
{
    // Polinomio es de la forma ax^3 + bx^2 + cx + d
    std::vector<double> resultado(2); // El resultado tendrá 2 coeficientes para un polinomio de grado 3

    resultado[1] = a;                    // El coeficiente de mayor grado se mantiene igual
    resultado[0] = b - r * resultado[1]; // Siguiente coeficiente

    // Invertir el vector del resultado
    std::reverse(resultado.begin(), resultado.end());

    return resultado;
}

std::vector<double> dividirPolinomios(double a, double b, double c, double r, double s) // division grado 4
{
    // Polinomio es de la forma ax^4 + bx^3 + cx^2 + dx + e
    // Divisor de la forma x^2 + bx + c

    std::vector<double> resultado(3); // El resultado tendrá 3 coeficientes para un polinomio de grado 4

    resultado[2] = a;                                       // El coeficiente de mayor grado se mantiene igual
    resultado[1] = b - r * resultado[2];                    // Siguiente coeficiente
    resultado[0] = c - r * resultado[1] - s * resultado[2]; // Siguiente coeficiente

    // Invertir el vector del resultado
    std::reverse(resultado.begin(), resultado.end());

    return resultado;
}

//--------------------- METODOS PARA ENCONTRAR RAICES DE POLINOMIOS -----------------------

void metodoBairstow(const std::vector<double> &coeficientes, int grado) // 1
{
    double ErrorR = 1000;
    double ErrorS = 1000;
    int iteraciones = 0;

    double r, s;
    std::cout << "Ingrese el valor inicial para r: ";
    std::cin >> r;
    std::cout << "Ingrese el valor inicial para s: ";
    std::cin >> s;
    std::cout << "\n";

    while (ErrorR > 1e-6 || ErrorS > 1e-6)
    {
        // double rAnterior = r; // Guardar el valor anterior de r
        // double sAnterior = s; // Guardar el valor anterior de s

        // Inicializar b4 dependiendo del grado
        double b4 = (grado == 3) ? 0 : coeficientes[0];

        double b3 = coeficientes[grado - 3] + b4 * r;
        double b2 = coeficientes[grado - 2] + r * b3 + s * b4;
        double b1 = coeficientes[grado - 1] + r * b2 + s * b3;
        double b0 = coeficientes[grado] + r * b1 + s * b2;

        double c4 = b4;
        double c3 = b3 + r * c4;
        double c2 = b2 + r * c3 + s * c4;
        double c1 = b1 + r * c2 + s * c3;

        double Ds_ = (c2 * c2) - (c1 * c3);
        double Dr = (-b1 * c2) - (-b0 * c3);
        double Ds = (c2 * -b0) - (c1 * -b1);

        double deltaR = Dr / Ds_;
        double deltaS = Ds / Ds_;

        r += deltaR;
        s += deltaS;

        ErrorR = std::abs(deltaR / r) * 100.0;
        ErrorS = std::abs(deltaS / s) * 100.0;

        iteraciones += 1;
        std::cout << "Iteracion actual: " << iteraciones << "\n";
        std::cout << "Valor de r: " << r << "\n";
        std::cout << "Valor de s: " << s << "\n";
        std::cout << "Porcentaje de error en r: " << ErrorR << "%\n";
        std::cout << "Porcentaje de error en s: " << ErrorS << "%\n";
        std::cout << "\n";
    }

    // Invertimos el signo de r y s
    r = -r;
    s = -s;

    auto raices = resolverx_2(1, r, s);

    std::cout << "Las raices del polinomio: ";
    imprimirPolinomio(grado, coeficientes);
    std::cout << " son:"
              << "\n";
    std::cout << "\n";
    std::cout << std::fixed << std::setprecision(30);

    std::cout << "Raiz 1 = " << raices.first << "\n";
    std::cout << "Raiz 2 = " << raices.second << "\n";

    if (grado == 3)
    {
        std::vector<double> resultado = dividirPolinomios(coeficientes[0], coeficientes[1], r);

        double Raiz3 = (-resultado[1]) / resultado[0];

        std::cout << "Raiz 3 = " << Raiz3 << "\n";
    }
    if (grado == 4)
    {
        std::vector<double> resultado = dividirPolinomios(coeficientes[0], coeficientes[1], coeficientes[2], r, s);

        double a = resultado[0];
        double b = resultado[1];
        double c = resultado[2];

        // std::cout << "a = " << a << "\n";
        // std::cout << "b = " << b << "\n";
        // std::cout << "c = " << c << "\n";

        auto raices = resolverx_2(a, b, c);

        std::cout << "Raiz 3 = " << raices.first << "\n";
        std::cout << "Raiz 4 = " << raices.second << "\n";
    }
}

void metodoCardano(const std::vector<double> &coeficientes, int grado) // 2
{
    if (grado == 3)
    {
        // double a = coeficientes[0];
        double b = coeficientes[1];
        double c = coeficientes[2];
        double d = coeficientes[3];

        double p = c - ((b * b) / 3);
        double q = (((2 * b * b * b) / 27) - ((b * c) / 3)) + d;

        double t_a = 1;
        double t_b = q;
        double t_c = -((p * p * p) / (27));

        auto raices_t = resolverx_2(t_a, t_b, t_c);

        double parteReal_t = raices_t.first.real();
        double parteImaginaria_t = raices_t.first.imag();

        double z = 0;
        double y = 0;
        double x = 0;

        if (parteImaginaria_t == 0)
        {

            z = cbrt(parteReal_t);
            y = z - (p / (3 * z));
        }
        else
        { // transformacion a coordenadas polares

            double z_modulo = sqrt((parteReal_t * parteReal_t) + (parteImaginaria_t * parteImaginaria_t));

            z = cbrt(z_modulo);

            double argumento_deg = atan2EnGrados(parteImaginaria_t, parteReal_t); // atan2 considera correctamente todos los cuadrantes

            double argumento = argumento_deg / 3;

            double z_real = z * cosenoEnGrados(argumento);
            double z_compleja = z * senoEnGrados(argumento);

            // declaramos el divisor real e imaginario para el conjugado

            double divisor_real = (3 * z_real);
            double divisor_complejo = (3 * z_compleja);

            double primer_termino = divisor_real * divisor_real;
            double segundo_termino = divisor_complejo * divisor_complejo;

            double divisor = primer_termino + segundo_termino; // obtenemos el divisor final

            double dividendo_real = p * divisor_real;
            double dividendo_complejo = p * -divisor_complejo;

            double valor_final_real = dividendo_real / divisor;
            double valor_final_complejo = dividendo_complejo / divisor;

            double y_real = z_real - (valor_final_real);
            double y_compleja = z_compleja - (valor_final_complejo);

            std::cout << "Las raices de t se han vuelto complejas y se han hecho las conversiones necesarias para continuar con la operacion\n";
            std::cout << "\n";
            std::cout << "Comprobando que el numero complejo se vuelve 0 (residuo de la resta): = " << y_compleja << "\n";
            std::cout << "\n";

            y = y_real;
        }

        x = y - (b / 3);

        std::vector<double> nuevoPolinomio = divisionSintetica(coeficientes[0], coeficientes[1], coeficientes[2], coeficientes[3], x);

        auto raices = resolverx_2(nuevoPolinomio[0], nuevoPolinomio[1], nuevoPolinomio[2]);

        std::cout << "Las raices del polinomio: ";
        imprimirPolinomio(grado, coeficientes);
        std::cout << " son:"
                  << "\n";
        std::cout << "\n";
        std::cout << std::fixed << std::setprecision(30);

        std::cout << "Raiz 1 = " << x << "\n";
        std::cout << "Raiz 2 = " << raices.first << "\n";
        std::cout << "Raiz 3 = " << raices.second << "\n";
    }

    if (grado == 4)
    {

        // obtenemos el polinomio dividiendo todo entre a (x^4 = 1)

        // double a_x = coeficientes[0] / coeficientes[0]; no se utiliza, ya que "a" siempre sera 1
        double b_x = coeficientes[1] / coeficientes[0];
        double c_x = coeficientes[2] / coeficientes[0];
        double d_x = coeficientes[3] / coeficientes[0];
        double e_x = coeficientes[4] / coeficientes[0];

        // primero comprobamos si es un tcp

        double comprobacion_ax = (((b_x * b_x) / 4) - c_x);
        double comprobacion_bx = -d_x; //-d
        double comprobacion_cx = -e_x; //-e

        if (esCuadradoPerfecto(comprobacion_ax, comprobacion_bx, comprobacion_cx))
        {
            std::cout << "La primera evaluacion arroja que tu polinomio es un cuadrado perfecto, se procedera de manera directa" << std::endl;
            std::cout << "\n";

            double derecho_x_1 = sqrt(comprobacion_ax);
            double derecho_x_0 = sqrt(comprobacion_cx) * obtenerSigno(comprobacion_bx);

            // std::cout << "lado derecho x^1: " << derecho_x_1 << "\n";
            // std::cout << "lado derecho x^0: " << derecho_x_0 << "\n";

            double izquierdo_x_2 = 1;
            double izquierdo_x_1 = b_x / 2;

            // std::cout << "lado izquierdo x^2: " << derecho_x_1 << "\n";
            // std::cout << "lado izquierdo x^1: " << derecho_x_0 << "\n";

            double nuevo_pol_x_2 = izquierdo_x_2;
            double nuevo_pol_x_1 = izquierdo_x_1 - derecho_x_1;
            double nuevo_pol_x_0 = -derecho_x_0;

            auto raices_finales = resolverx_2(nuevo_pol_x_2, nuevo_pol_x_1, nuevo_pol_x_0);

            std::cout << "Las raices del polinomio: ";
            imprimirPolinomio(grado, coeficientes);
            std::cout << " son:"
                      << "\n";
            std::cout << "\n";
            std::cout << std::fixed << std::setprecision(30);

            std::cout << "Raiz 1 = " << raices_finales.first << "\n";
            std::cout << "Raiz 2 = " << raices_finales.second << "\n";

            std::vector<double> nuevo_polinomio = dividirPolinomios(coeficientes[0], coeficientes[1], coeficientes[2], nuevo_pol_x_1, nuevo_pol_x_0);

            double a = nuevo_polinomio[0];
            double b = nuevo_polinomio[1];
            double c = nuevo_polinomio[2];

            auto raices = resolverx_2(a, b, c);

            std::cout << "Raiz 3 = " << raices.first << "\n";
            std::cout << "Raiz 4 = " << raices.second << "\n";
        }
        else
        {
            std::cout << "La primera evaluacion arroja que tu polinomio no es un cuadrado perfecto, se convertira a uno en este momento" << std::endl;
            std::cout << "\n";

            std::vector<double> PolinomioTartaglia(3);

            // sumas a los coeficientes

            double b_y = c_x * (-1);
            double c_y = ((4 * e_x) - (b_x * d_x)) * (-1);
            double d_y = ((b_x * b_x * e_x) - (4 * c_x * e_x) + (d_x * d_x)) * (-1);

            double p = c_y - ((b_y * b_y) / 3);
            double q = (((2 * b_y * b_y * b_y) / 27) - ((b_y * c_y) / 3)) + d_y;

            double t_a = 1;
            double t_b = q;
            double t_c = -((p * p * p) / (27));

            auto raices_t = resolverx_2(t_a, t_b, t_c);

            double parteReal_t = raices_t.first.real();
            double parteImaginaria_t = raices_t.first.imag();

            double z = 0;
            double y = 0;
            double y_prima = 0;

            if (parteImaginaria_t == 0)
            {
                // solo necesitamos una raiz

                z = cbrt(parteReal_t);
                y_prima = z - (p / (3 * z));
            }
            else // coordenadas polares cuando t es complejo
            {
                double z_modulo = sqrt((parteReal_t * parteReal_t) + (parteImaginaria_t * parteImaginaria_t));

                z = cbrt(z_modulo);

                double argumento_deg = atan2EnGrados(parteImaginaria_t, parteReal_t); // atan2 considera correctamente todos los cuadrantes

                double argumento = argumento_deg / 3;

                double z_real = z * cosenoEnGrados(argumento);
                double z_compleja = z * senoEnGrados(argumento);

                // declaramos el divisor real e imaginario para el conjugado

                double divisor_real = (3 * z_real);
                double divisor_complejo = (3 * z_compleja);

                double primer_termino = divisor_real * divisor_real;
                double segundo_termino = divisor_complejo * divisor_complejo;

                double divisor = primer_termino + segundo_termino; // obtenemos el divisor final

                double dividendo_real = p * divisor_real;
                double dividendo_complejo = p * -divisor_complejo;

                double valor_final_real = dividendo_real / divisor;
                double valor_final_complejo = dividendo_complejo / divisor;

                double y_real = z_real - (valor_final_real);
                double y_compleja = z_compleja - (valor_final_complejo);

                std::cout << "Las raices de t se han vuelto complejas y se han hecho las conversiones necesarias para continuar con la operacion\n";
                std::cout << "\n";
                std::cout << "Comprobando que el numero complejo se vuelve 0 (residuo de la resta): = " << y_compleja << "\n";
                std::cout << "\n";

                y_prima = y_real;
            }

            // std::cout << "la y prima resultante= " << y_prima << "\n";
            // std::cout << "\n";

            y = y_prima - (b_y / 3);

            // auto fraccion = decimalAFraccion(y);
            // std::cout << fraccion.numerator() << "/" << fraccion.denominator() << std::endl;

            // lado izquierdo de la ecuacion

            double coef_x_2i = 1;
            double coef_x_1i = b_x / 2;
            double coef_x_0i = y / 2;

            // lado derecho de la ecuacion

            double coef_x_2d = ((b_x * b_x) / 4) - c_x + y;
            double coef_x_1d = ((b_x * y) / 2) - d_x;
            double coef_x_0d = ((y * y) / 4) - e_x;

            // verificar si es tcp

            if (esCuadradoPerfecto(coef_x_2d, coef_x_1d, coef_x_0d))
            {
                std::cout << "El polinomio se ha convertido a cuadrado perfecto exitosamente!" << std::endl;
                std::cout << "\n";

                // nuevos coeficientes para el lado derecho

                // double tolerancia = 1e-6;

                double nuevo_x_1 = sqrt(coef_x_2d);
                double nuevo_x_0 = sqrt(coef_x_0d) * obtenerSigno(coef_x_1d);
                ;

                double nuevo_pol_x_2 = coef_x_2i;
                double nuevo_pol_x_1 = (coef_x_1i) - (nuevo_x_1);
                double nuevo_pol_x_0 = (coef_x_0i) - (nuevo_x_0);

                // std::cout << "b_x " << b_x << "\n";
                // std::cout << "c_x " << c_x << "\n";
                // std::cout << "b_x " << b_x << "\n";

                // std::cout << "coeficiente de x cuadrada para el nuevo polinomio " << coef_x_2d << "\n";
                // std::cout << "coeficiente de x1 para el nuevo polinomio " << coef_x_1d << "\n";
                // std::cout << "coeficiente de x0 para el nuevo polinomio " << coef_x_0d << "\n";

                auto raices_finales = resolverx_2(nuevo_pol_x_2, nuevo_pol_x_1, nuevo_pol_x_0);

                std::cout << "Las raices del polinomio: ";
                imprimirPolinomio(grado, coeficientes);
                std::cout << " son:"
                          << "\n";
                std::cout << "\n";
                std::cout << std::fixed << std::setprecision(30);

                std::cout << "Raiz 1 = " << raices_finales.first << "\n";
                std::cout << "Raiz 2 = " << raices_finales.second << "\n";

                std::vector<double> nuevo_polinomio = dividirPolinomios(coeficientes[0], coeficientes[1], coeficientes[2], nuevo_pol_x_1, nuevo_pol_x_0);

                double a = nuevo_polinomio[0];
                double b = nuevo_polinomio[1];
                double c = nuevo_polinomio[2];

                auto raices = resolverx_2(a, b, c);

                std::cout << "Raiz 3 = " << raices.first << "\n";
                std::cout << "Raiz 4 = " << raices.second << "\n";
            }
            else
            {
                std::cout << "El polinomio no es un cuadrado perfecto. Verifica los calculos internos" << std::endl;
            }
        }
    }
}

void metodoBiseccion(const std::vector<double> &coeficientes, int grado) // 3
{
    if (grado == 4)
    {

        double intervaloA, intervaloB;
        std::cout << "Ingrese el valor del intervalo [a]: ";
        std::cin >> intervaloA;
        std::cout << "Ingrese el valor del intervalo [b]: ";
        std::cin >> intervaloB;
        std::cout << "\n";

        double f_intervaloA = (coeficientes[0] * intervaloA * intervaloA * intervaloA * intervaloA) + (coeficientes[1] * intervaloA * intervaloA * intervaloA) +
                              (coeficientes[2] * intervaloA * intervaloA) + (coeficientes[3] * intervaloA) + coeficientes[4];
        double f_intervaloB = (coeficientes[0] * intervaloA * intervaloA * intervaloA * intervaloA) + (coeficientes[1] * intervaloA * intervaloA * intervaloA) +
                              (coeficientes[2] * intervaloA * intervaloA) + (coeficientes[3] * intervaloA) + coeficientes[4];

        if (f_intervaloA * f_intervaloB >= 0)
        {
            std::cout << "ATENCION: Puede que el intervalo que introdujiste no converga correctamente!\n";
            std::cout << "\n";
        }

        double Error = 1000;
        int iteraciones = 0;

        double pn{}, f_pn;

        while (Error > 1e-15)
        {
            double pn_anterior = pn;
            pn = (intervaloA + intervaloB) / 2; // punto medio
            f_pn = (coeficientes[0] * pn * pn * pn * pn) + (coeficientes[1] * pn * pn * pn) +
                   (coeficientes[2] * pn * pn) + (coeficientes[3] * pn) + coeficientes[4];

            if (f_pn * f_intervaloA < 0) // establecemos condicion para saber de que lado colocar pn
            {
                intervaloB = pn;
            }
            else
            {
                intervaloA = pn;
            }

            Error = std::abs((pn - pn_anterior) / pn);

            iteraciones += 1;

            std::cout << "Iteracion actual: " << iteraciones << "\n";
            std::cout << "Raiz aproximada: " << pn << "\n";
            std::cout << "Valor de f(p): " << f_pn << "\n";
            std::cout << "Porcentaje de error: " << Error * 100.0 << "%\n";
            std::cout << "\n";

            if (iteraciones > 100 && std::abs(Error) > 1e-6)
            {
                std::cerr << "Se han realizado: " << iteraciones << " iteraciones, y el metodo no converge correctamente. Revisa tu intervalo." << std::endl;
                std::cout << "\n";
                return;
            }
        }

        double comprobacion = evaluarPolinomio(coeficientes, pn);
        // std::cout << "comprobacion" << comprobacion << "\n";

        if (std::abs(comprobacion) < 1e-6)
        {
            std::vector<double> nuevoPolinomio = divisionSintetica(coeficientes[0], coeficientes[1], coeficientes[2], coeficientes[3], coeficientes[4], pn);

            std::cout << "Se procedera con Tartaglia Cardano para obtener raices del nuevo polinomio de grado 3: \n";
            std::cout << "\n";

            std::cout << "Las raices del polinomio: ";
            imprimirPolinomio(grado, coeficientes);
            std::cout << " son:"
                      << "\n";
            std::cout << "\n";
            std::cout << std::fixed << std::setprecision(30);

            metodoCardano(nuevoPolinomio, 3);
            std::cout << "Raiz 4 = " << pn << "\n";
            std::cout << "\n";
            return;
        }
        else
        {
            system("cls");
            std::cout << "Comprobando la raiz, se ha calculado que el valor no ha podido converger, revisa tus valores e intentalo de nuevo!\n";
            std::cout << "\n";
            std::cout << "Raiz encontrada: " << pn << "\n";
            std::cout << "\n";
        }
    }

    if (grado == 3)
    {

        double intervaloA, intervaloB;
        std::cout << "Ingrese el valor del intervalo [a]: ";
        std::cin >> intervaloA;
        std::cout << "Ingrese el valor del intervalo [b]: ";
        std::cin >> intervaloB;
        std::cout << "\n";

        double f_intervaloA = (coeficientes[0] * intervaloA * intervaloA * intervaloA) +
                              (coeficientes[1] * intervaloA * intervaloA) + (coeficientes[2] * intervaloA) + coeficientes[3];
        double f_intervaloB = (coeficientes[0] * intervaloA * intervaloA * intervaloA) +
                              (coeficientes[1] * intervaloA * intervaloA) + (coeficientes[2] * intervaloA) + coeficientes[3];

        if (f_intervaloA * f_intervaloB >= 0)
        {
            std::cout << "ATENCION: Puede que el intervalo que introdujiste no converga correctamente!\n";
            std::cout << "\n";
        }

        double Error = 1000;
        int iteraciones = 0;

        double pn = 0.0, f_pn;

        while (Error > 1e-15)
        {
            double pn_anterior = pn;
            pn = (intervaloA + intervaloB) / 2; // punto medio
            f_pn = (coeficientes[0] * pn * pn * pn) +
                   (coeficientes[1] * pn * pn) + (coeficientes[2] * pn) + coeficientes[3];

            if (f_pn * f_intervaloA < 0) // establecemos condicion para saber de que lado colocar pn
            {
                intervaloB = pn;
            }
            else
            {
                intervaloA = pn;
            }

            Error = std::abs(pn - pn_anterior) / std::abs(pn);

            iteraciones += 1;

            std::cout << "Iteracion actual: " << iteraciones << "\n";
            std::cout << "Raiz aproximada: " << pn << "\n";
            std::cout << "Valor de f(p): " << f_pn << "\n";
            std::cout << "Porcentaje de error: " << Error * 100.0 << "%\n";
            std::cout << "\n";
        }

        double comprobacion = evaluarPolinomio(coeficientes, pn);
        // std::cout << "comprobacion" << comprobacion << "\n";

        if (std::abs(comprobacion) < 1e-6)
        {
            std::cout << "Las raices del polinomio: ";
            imprimirPolinomio(grado, coeficientes);
            std::cout << " son:"
                      << "\n";
            std::cout << "\n";
            std::cout << std::fixed << std::setprecision(30);

            std::vector<double> nuevoPolinomio = divisionSintetica(coeficientes[0], coeficientes[1], coeficientes[2], coeficientes[3], coeficientes[4], pn);

            auto raices = resolverx_2(nuevoPolinomio[0], nuevoPolinomio[1], nuevoPolinomio[2]);

            std::cout << "Raiz 1 = " << pn << "\n";
            std::cout << "Raiz 2 = " << raices.first << "\n";
            std::cout << "Raiz 3 = " << raices.second << "\n";
            std::cout << "\n";
        }
        else
        {
            system("cls");
            std::cout << "Comprobando la raiz, se ha calculado que el valor no ha podido converger, revisa tus valores e intentalo de nuevo!\n";
            std::cout << "\n";
            std::cout << "Raiz encontrada: " << pn << "\n";
            std::cout << "\n";
        }
    }
}

void metodoMuller(const std::vector<double> &coeficientes, int grado) // 4
{
    double error = 1000;
    int iteraciones = 0;

    double x0, x1, x2, x3;
    double fx0, fx1, fx2;

    double h0, h1;
    double delta0, delta1;

    double a, b, c;

    std::cout << "Ingrese la primera condicion inicial (x0): ";
    std::cin >> x0;
    std::cout << "Ingrese la segunda condicion inicial (x1): ";
    std::cin >> x1;
    std::cout << "Ingrese la tercera condicion inicial (x2): ";
    std::cin >> x2;
    std::cout << "\n";

    while (error > 1e-6)
    {
        // evaluar el polinomio

        fx0 = evaluarPolinomio(coeficientes, x0);
        fx1 = evaluarPolinomio(coeficientes, x1);
        fx2 = evaluarPolinomio(coeficientes, x2);

        // valores de las iteraciones

        h0 = x1 - x0;
        h1 = x2 - x1;

        delta0 = (fx1 - fx0) / h0;
        delta1 = (fx2 - fx1) / h1;

        a = (delta1 - delta0) / (h1 + h0);
        b = (a * h1) + delta1;
        c = fx2;

        x3 = x2 - ((2 * c) / (b + ((obtenerSigno(b)) * (sqrt((b * b) - 4 * a * c)))));

        if (std::isnan(x3))
        {
            std::cout << "Ha ocurrido un error con los calculos, verifica las condiciones iniciales establecidas y vuelve a intentarlo." << std::endl;
            std::cout << "\n";
            break;
        }

        error = std::abs((x3 - x2) / x3) * 100.0;

        iteraciones += 1;

        std::cout << "Iteracion actual: " << iteraciones << "\n";
        std::cout << "Valor de a: " << a << "\n";
        std::cout << "Valor de b: " << b << "\n";
        std::cout << "Valor de c: " << c << "\n";
        std::cout << "Valor de x (x3): " << x3 << "\n";
        std::cout << "Porcentaje de error: " << error << "%\n";
        std::cout << "\n";

        x0 = x1;
        x1 = x2;
        x2 = x3;
    }

    if (grado == 3 && !std::isnan(x3))
    {
        std::vector<double> nuevoPolinomio = divisionSintetica(coeficientes[0], coeficientes[1], coeficientes[2], coeficientes[3], coeficientes[4], x3);

        auto raices = resolverx_2(nuevoPolinomio[0], nuevoPolinomio[1], nuevoPolinomio[2]);

        std::cout << "Las raices del polinomio: ";
        imprimirPolinomio(grado, coeficientes);
        std::cout << " son:"
                  << "\n";
        std::cout << "\n";
        std::cout << std::fixed << std::setprecision(30);

        std::cout << "Raiz 1 = " << x3 << "\n";
        std::cout << "Raiz 2 = " << raices.first << "\n";
        std::cout << "Raiz 3 = " << raices.second << "\n";
    }
    if (grado == 4 && !std::isnan(x3))
    {
        std::vector<double> nuevoPolinomio = divisionSintetica(coeficientes[0], coeficientes[1], coeficientes[2], coeficientes[3], coeficientes[4], x3);

        std::cout << "Se procedera con Tartaglia Cardano para obtener raices del nuevo polinomio de grado 3: \n";
        std::cout << "\n";

        std::cout << "Las raices del polinomio: ";
        imprimirPolinomio(grado, coeficientes);
        std::cout << " son:"
                  << "\n";
        std::cout << "\n";
        std::cout << std::fixed << std::setprecision(30);

        metodoCardano(nuevoPolinomio, 3);
        std::cout << "Raiz 4 = " << x3 << "\n";
    }
}

void metodoNewton(const std::vector<double> &coeficientes, int grado) // 5
{
    double intervaloA, intervaloB;
    std::cout << "Ingrese el valor del intervalo [a]: ";
    std::cin >> intervaloA;
    std::cout << "Ingrese el valor del intervalo [b]: ";
    std::cin >> intervaloB;
    std::cout << "\n";

    double x = (intervaloA + intervaloB) / 2; // calculamos el punto medio

    // primer evaluacion:

    double fx = evaluarPolinomio(coeficientes, x);
    double dfx = evaluarDerivada(coeficientes, x);

    double error = 1000;
    int iteraciones = 0;

    while (std::abs(fx) > 1e-14 && std::abs(error) > 1e-14) // cuando el porcentaje de error disminuya lo suficiente
    {
        double x_anterior = x;

        if (dfx == 0.0)
        {
            std::cerr << "La derivada es cero en el punto medio obtenido: " << x << ", no se pudo calcular el resultado!" << std::endl;
            return;
        }

        x = x - (fx / dfx); // fórmula de Newton-Raphson

        fx = evaluarPolinomio(coeficientes, x); // Recalcular fx
        dfx = evaluarDerivada(coeficientes, x); // Recalcular dfx

        error = (x - x_anterior) / x * 100;

        iteraciones += 1;

        std::cout << "Iteracion actual: " << iteraciones << "\n";
        std::cout << "Valor de x actual: " << x << "\n";
        std::cout << "Porcentaje de error: " << error << "%\n";
        std::cout << "\n";

        if (iteraciones > 100 && std::abs(error) > 1e-6)
        {
            std::cerr << "Se han realizado: " << iteraciones << " iteraciones, y el metodo no converge correctamente. Revisa tu intervalo." << std::endl;
            std::cout << "\n";
            return;
        }
    }

    double comprobacion = evaluarPolinomio(coeficientes, x);
    // std::cout << "comprobacion" << comprobacion << "\n";

    if (std::abs(comprobacion) < 1e-6)
    {
        if (grado == 4)
        {
            std::vector<double> nuevoPolinomio = divisionSintetica(coeficientes[0], coeficientes[1], coeficientes[2], coeficientes[3], coeficientes[4], x);

            std::cout << "Se procedera con Tartaglia Cardano para obtener raices del nuevo polinomio de grado 3: \n";
            std::cout << "\n";

            std::cout << "Las raices del polinomio: ";
            imprimirPolinomio(grado, coeficientes);
            std::cout << " son:"
                      << "\n";
            std::cout << "\n";
            std::cout << std::fixed << std::setprecision(30);

            metodoCardano(nuevoPolinomio, 3);
            std::cout << "Raiz 4 = " << x << "\n";
        }
        if (grado == 3)
        {
            std::vector<double> nuevoPolinomio = divisionSintetica(coeficientes[0], coeficientes[1], coeficientes[2], coeficientes[3], coeficientes[4], x);

            auto raices = resolverx_2(nuevoPolinomio[0], nuevoPolinomio[1], nuevoPolinomio[2]);

            std::cout << "Las raices del polinomio: ";
            imprimirPolinomio(grado, coeficientes);
            std::cout << " son:"
                      << "\n";
            std::cout << "\n";
            std::cout << std::fixed << std::setprecision(30);

            std::cout << "Raiz 1 = " << x << "\n";
            std::cout << "Raiz 2 = " << raices.first << "\n";
            std::cout << "Raiz 3 = " << raices.second << "\n";
        }
    }
    else
    {
        system("cls");
        std::cout << "Comprobando la raiz, se ha calculado que el valor no ha podido converger, revisa tus valores e intentalo de nuevo!\n";
        std::cout << "\n";
        std::cout << "Raiz encontrada: " << x << "\n";
        std::cout << "\n";
    }
}

//------------------- METODOS PARA ENCONTRAR RAICES DE SISTEMAS DE ECUACIONES --------------

void metodoGaussSeidel(const std::vector<double> &CoeficientesX, const std::vector<double> &CoeficientesY, const std::vector<double> &CoeficientesZ, const std::vector<double> &terminosIndependientes)
{
    double x, y, z;
    double errorX, errorY, errorZ;

    // condiciones iniciales
    std::cout << "Ingrese la primera condicion inicial (x): ";
    std::cin >> x;
    std::cout << "Ingrese la segunda condicion inicial (y): ";
    std::cin >> y;
    std::cout << "Ingrese la tercera condicion inicial (z): ";
    std::cin >> z;
    std::cout << "\n";

    double precision = 1e-9;
    int iteraciones = 0;
    errorX = errorY = errorZ = 1000;

    while (std::abs(errorX) > precision || std::abs(errorY) > precision || std::abs(errorZ) > precision)
    {
        // valores anteriores (% error)
        double x_anterior = x;
        double y_anterior = y;
        double z_anterior = z;

        x = ((-CoeficientesY[0] * y) + (-CoeficientesZ[0] * z) + terminosIndependientes[0]) / CoeficientesX[0];
        y = ((-CoeficientesX[1] * x) + (-CoeficientesZ[1] * z) + terminosIndependientes[1]) / CoeficientesY[1];
        z = ((-CoeficientesX[2] * x) + (-CoeficientesY[2] * y) + terminosIndependientes[2]) / CoeficientesZ[2];

        errorX = (x - x_anterior) / x * 100;
        errorY = (y - y_anterior) / y * 100;
        errorZ = (z - z_anterior) / z * 100;

        iteraciones += 1;

        std::cout << "Iteracion actual: " << iteraciones << "\n";
        std::cout << "Valor de x actual: " << x << "\n";
        std::cout << "Valor de y actual: " << y << "\n";
        std::cout << "Valor de z actual: " << z << "\n";
        std::cout << "Porcentaje de error en x: " << errorX << "%\n";
        std::cout << "Porcentaje de error en y: " << errorY << "%\n";
        std::cout << "Porcentaje de error en z: " << errorZ << "%\n";
        std::cout << "\n";

        if (iteraciones > 100 && std::abs(errorX) > 1e-6 && std::abs(errorY) > 1e-6 && std::abs(errorZ) > 1e-6)
        {
            system("cls");
            std::cerr << "Se han realizado: " << iteraciones << " iteraciones, y el metodo no converge correctamente. Revisa los coeficientes o el orden de las ecuaciones y vuelve a intentarlo." << std::endl;
            std::cout << "\n";
            std::cout << "Valor de x encontrado: " << x << "\n";
            std::cout << "Valor de y encontrado: " << y << "\n";
            std::cout << "Valor de z encontrado: " << z << "\n";
            std::cout << "\n";
            return;
        }
    }

    std::cout << "El metodo logro converger despues de " << iteraciones << " iteraciones." << std::endl;
    std::cout << "Solucion aproximada: x = " << x << std::endl;
    std::cout << "Solucion aproximada: y = " << y << std::endl;
    std::cout << "Solucion aproximada: z = " << z << std::endl;
    std::cout << "\n";
}

void metodoJacobi(const std::vector<double> &CoeficientesX, const std::vector<double> &CoeficientesY, const std::vector<double> &CoeficientesZ, const std::vector<double> &terminosIndependientes)
{
    double x, y, z;
    double errorX, errorY, errorZ;

    // condiciones iniciales
    std::cout << "Ingrese la primera condicion inicial (x): ";
    std::cin >> x;
    std::cout << "Ingrese la segunda condicion inicial (y): ";
    std::cin >> y;
    std::cout << "Ingrese la tercera condicion inicial (z): ";
    std::cin >> z;
    std::cout << "\n";

    double precision = 1e-9;
    int iteraciones = 0;
    errorX = errorY = errorZ = 1000;

    while (std::abs(errorX) > precision || std::abs(errorY) > precision || std::abs(errorZ) > precision)
    {
        // valores anteriores (% error)
        double x_anterior = x;
        double y_anterior = y;
        double z_anterior = z;

        x = ((-CoeficientesY[0] * y_anterior) + (-CoeficientesZ[0] * z_anterior) + terminosIndependientes[0]) / CoeficientesX[0];
        y = ((-CoeficientesX[1] * x_anterior) + (-CoeficientesZ[1] * z_anterior) + terminosIndependientes[1]) / CoeficientesY[1];
        z = ((-CoeficientesX[2] * x_anterior) + (-CoeficientesY[2] * y_anterior) + terminosIndependientes[2]) / CoeficientesZ[2];

        errorX = (x - x_anterior) / x * 100;
        errorY = (y - y_anterior) / y * 100;
        errorZ = (z - z_anterior) / z * 100;

        iteraciones += 1;

        std::cout << "Iteracion actual: " << iteraciones << "\n";
        std::cout << "Valor de x actual: " << x << "\n";
        std::cout << "Valor de y actual: " << y << "\n";
        std::cout << "Valor de z actual: " << z << "\n";
        std::cout << "Porcentaje de error en x: " << errorX << "%\n";
        std::cout << "Porcentaje de error en y: " << errorY << "%\n";
        std::cout << "Porcentaje de error en z: " << errorZ << "%\n";
        std::cout << "\n";

        if (iteraciones >= 100 && std::abs(errorX) > 1e-6 && std::abs(errorY) > 1e-6 && std::abs(errorZ) > 1e-6)
        {
            system("cls");
            std::cerr << "Se han realizado" << iteraciones << "iteraciones, y el metodo no converge correctamente. Revisa los coeficientes o el orden de las ecuaciones y vuelve a intentarlo." << std::endl;
            std::cout << "\n";
            std::cout << "Valor de x encontrado: " << x << "\n";
            std::cout << "Valor de y encontrado: " << y << "\n";
            std::cout << "Valor de z encontrado: " << z << "\n";
            std::cout << "\n";
            return;
        }
    }

    std::cout << "El metodo logro converger despues de " << iteraciones << " iteraciones." << std::endl;
    std::cout << "Solucion aproximada: x = " << x << std::endl;
    std::cout << "Solucion aproximada: y = " << y << std::endl;
    std::cout << "Solucion aproximada: z = " << z << std::endl;
    std::cout << "\n";
}

//--------------------------------------- MAIN ---------------------------------------------

int main()
{
    int opcion;
    char continuar;

    do
    {
        std::cout << "Calculadora de Metodos Numericos" << std::endl;
        std::cout << "\n";
        std::cout << "Que deseas calcular?" << std::endl;
        std::cout << "\n";
        std::cout << "1) Polinomios" << std::endl;
        std::cout << "2) Sistema de ecuaciones (3 incognitas)" << std::endl;
        std::cout << "3) Salir" << std::endl;
        std::cout << "\n";
        std::cout << "Selecciona una opcion: ";
        std::cin >> opcion;
        system("cls");

        switch (opcion)
        {
        case 1:
            int grado;
            std::cout << "Ingrese el grado del polinomio: ";
            std::cin >> grado;
            std::cout << "\n";

            if (grado > 4 || grado < 1)
            {
                std::cout << "La calculadora solo maneja hasta polinomios de grado 4, si deseas calcular las raices de polinomios de un grado superior, es recomendable utilizar mi calculadora de Python :)\n";
                std::cout << "\n";
            }
            else
            {
                std::vector<double> coeficientes(grado + 1);
                for (int i = grado; i >= 0; i--)
                {
                    std::cout << "Ingrese el coeficiente de x^" << i << ": ";
                    std::cin >> coeficientes[grado - i];
                }

                system("cls");

                if (grado == 2)
                {
                    std::cout << std::fixed << std::setprecision(30);
                    std::cout << "Tu polinomio es de segundo grado, por formula general se ha obtenido: \n";
                    std::cout << "\n";

                    auto raices = resolverx_2(coeficientes[0], coeficientes[1], coeficientes[2]);

                    std::cout << "Raiz 1 = " << raices.first << "\n";
                    std::cout << "Raiz 2 = " << raices.second << "\n";
                }

                else if (grado == 1)
                {
                    std::cout << std::fixed << std::setprecision(30);
                    double raiz = resolverx_1(coeficientes[0], coeficientes[1]);

                    if (coeficientes[0] == 0)
                    {
                        std::cout << "No se pudo realizar el despeje, ya que el coeficiente de X no puede estar vacio, vuelve a ingresar los datos\n";
                        std::cout << "\n";
                        std::cout << "El numero que ingresaste es: " << coeficientes[1] << "\n";
                    }
                    else
                    {
                        std::cout << "La raiz encontrada es: = " << raiz << "\n";
                    }
                }

                else if (grado > 2 && grado <= 4)
                {
                    int opcion;

                    std::cout << "Bienvenido a la Calculadora de Polinomios\n";
                    std::cout << "\n";
                    std::cout << "El polinomio que has ingresado es:\n";
                    std::cout << "\n";

                    imprimirPolinomio(grado, coeficientes);

                    std::cout << std::fixed << std::setprecision(30);

                    std::cout << "\n";
                    std::cout << "\n";

                    std::cout << "1) Metodo de Bairstow\n";
                    std::cout << "2) Metodo de Tartaglia Cardano\n";
                    std::cout << "3) Metodo de Biseccion\n";
                    std::cout << "4) Metodo de Muller\n";
                    std::cout << "5) Metodo de Newton-Raphson\n";

                    std::cout << "\n";
                    std::cout << "Seleccione una opcion: ";
                    std::cin >> opcion;
                    std::cout << "\n";

                    switch (opcion)
                    {
                    case 1:
                        system("cls");
                        metodoBairstow(coeficientes, grado);
                        break;
                    case 2:
                        system("cls");
                        metodoCardano(coeficientes, grado);
                        break;
                    case 3:
                        system("cls");
                        metodoBiseccion(coeficientes, grado);
                        break;
                    case 4:
                        system("cls");
                        metodoMuller(coeficientes, grado);
                        break;
                    case 5:
                        system("cls");
                        metodoNewton(coeficientes, grado);
                        break;
                    default:
                        system("cls");
                        std::cout << "Opcion no valida.\n";
                        break;
                    }
                }
            }

            break;

        case 2:
        {
            std::cout << std::fixed << std::setprecision(2);
            std::vector<double> CoeficientesX(3), CoeficientesY(3), CoeficientesZ(3), terminosIndependientes(3);

            for (int i = 0; i < 3; i++)
            {
                std::cout << "ECUACION " << i + 1 << "\n";
                std::cout << "\n";
                std::cout << "Elija el valor de x: ";
                std::cin >> CoeficientesX[i];
                std::cout << "Elija el valor de y: ";
                std::cin >> CoeficientesY[i];
                std::cout << "Elija el valor de z: ";
                std::cin >> CoeficientesZ[i];
                std::cout << "Termino independiente: ";
                std::cin >> terminosIndependientes[i];
                std::cout << "\n";
            }

            system("cls");
            std::cout << "El sistema de ecuaciones que ingresaste fue:\n";
            std::cout << "\n";

            for (int i = 0; i < 3; i++)
            {
                std::cout << (CoeficientesX[i] >= 0 ? "+" : "") << CoeficientesX[i] << "x ";
                std::cout << (CoeficientesY[i] >= 0 ? "+" : "") << CoeficientesY[i] << "y ";
                std::cout << (CoeficientesZ[i] >= 0 ? "+" : "") << CoeficientesZ[i] << "z ";
                std::cout << "= " << terminosIndependientes[i] << "\n";
            }

            std::cout << std::fixed << std::setprecision(30);
            std::cout << "\n";
            std::cout << "\n";
            std::cout << "1) Metodo de Gauss-Seidel\n";
            std::cout << "2) Metodo de Jacobi\n";
            std::cout << "\n";
            std::cout << "Seleccione una opcion: ";
            std::cin >> opcion;
            std::cout << "\n";

            switch (opcion)
            {
            case 1:
                system("cls");
                metodoGaussSeidel(CoeficientesX, CoeficientesY, CoeficientesZ, terminosIndependientes);
                break;
            case 2:
                system("cls");
                metodoJacobi(CoeficientesX, CoeficientesY, CoeficientesZ, terminosIndependientes);
                break;
            default:
                system("cls");
                std::cout << "Opcion no valida.\n";
                break;
            }
            break; // Salir del switch exterior
        }
        case 3:
            std::cout << "Gracias por usar la calculadora! Elaborada por: Luna Cervantes Bernardo." << std::endl;
            return 0;
        default:
            std::cout << "Opcion no valida. Intentalo de nuevo!" << std::endl;
            std::cout << "\n";
            return 0;
        }

        std::cout << "Deseas realizar otra operacion? (s/n): ";
        std::cin >> continuar;
        system("cls");

    } while (continuar == 'S' || continuar == 's');

    std::cout << "Gracias por usar la calculadora! Elaborada por: Luna Cervantes Bernardo." << std::endl;
    return 0;
}

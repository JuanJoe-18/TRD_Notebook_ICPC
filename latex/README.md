# Compilar el notebook

El PDF se genera con XeLaTeX (usado por el libro: `latex/plantilla.tex`).

```sh
cd latex
xelatex plantilla.tex
xelatex plantilla.tex
```

Salida: `latex/plantilla.pdf`.

## Notas

- Hay que compilar **desde `latex/`** porque los `\lstinputlisting{../...}` son rutas relativas al directorio de trabajo.
- Se compila **dos veces**: el índice de contenidos resuelve las páginas en la segunda pasada.
- Requiere MiKTeX/TeXLive con `xelatex` y la fuente **Space Mono** (Google Fonts). Si no la tienes: instálala (Windows: `C:\Windows\Fonts` o por usuario en `%LOCALAPPDATA%\Microsoft\Windows\Fonts`). En sistemas sin ella se usa el fallback de `pdflatex` (Inconsolata).
- Alternativa con pdflatex: `pdflatex plantilla.tex` (misma salida, también dos pasadas).
- El encabezado ICPC (universidad/equipo y número de página) se configura en `plantilla.tex` (macros `\fancyhead`).
- Actualiza el PDF tras modificar cualquier `.hpp`: solo re-ejecuta los dos comandos.
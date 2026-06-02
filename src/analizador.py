import argparse
import sys


def leer_fasta(ruta):
    """
    Lee un archivo FASTA y extrae todas las secuencias.

    Args:
        ruta (str): Ruta del archivo FASTA.

    Returns:
        list: Lista de tuplas (encabezado, secuencia).
    """
    archivo = open(ruta, "r")
    lineas = archivo.readlines()
    archivo.close()

    secuencias = []
    encabezado = None
    secuencia = ""

    for linea in lineas:
        linea = linea.strip()
        if linea.startswith(">"):
            if encabezado is not None:
                secuencias.append((encabezado, secuencia))
            encabezado = linea
            secuencia = ""
        else:
            secuencia += linea

    if encabezado is not None:
        secuencias.append((encabezado, secuencia))

    return secuencias


def calcular_gc(secuencia):
    """
    Calcula el contenido GC de una secuencia de ADN.

    Args:
        secuencia (str): Secuencia de ADN.

    Returns:
        float: Proporción de bases G y C entre 0 y 1.
    """
    secuencia = secuencia.upper()
    g = secuencia.count("G")
    c = secuencia.count("C")
    gc = g + c

    if len(secuencia) == 0:
        return 0.0

    porcentaje = gc / len(secuencia)
    return porcentaje


def calcular_estadisticas(secuencias):
    """
    Calcula estadísticas para todas las secuencias.

    Args:
        secuencias (list): Lista de tuplas (encabezado, secuencia).

    Returns:
        list: Lista de diccionarios con estadísticas de cada secuencia.
    """
    estadisticas = []

    for encabezado, secuencia in secuencias:
        longitud = len(secuencia)
        gc = calcular_gc(secuencia)
        stats = {"encabezado": encabezado, "longitud": longitud, "contenido_gc": gc}
        estadisticas.append(stats)

    return estadisticas


def pasa_filtros(stats, min_len, max_len, min_gc, max_gc):
    """
    Verifica si una secuencia cumple todos los filtros.

    Args:
        stats (dict): Diccionario con estadísticas de la secuencia.
        min_len (int): Longitud mínima permitida.
        max_len (int): Longitud máxima permitida.
        min_gc (float): Contenido GC mínimo permitido.
        max_gc (float): Contenido GC máximo permitido.

    Returns:
        bool: True si cumple todos los filtros, False en caso contrario.
    """
    if stats["longitud"] < min_len:
        return False

    if stats["longitud"] > max_len:
        return False

    if stats["contenido_gc"] < min_gc:
        return False

    if stats["contenido_gc"] > max_gc:
        return False

    return True


def escribir_resultados(estadisticas_filtradas, ruta):
    """
    Escribe los resultados filtrados en un archivo TSV.

    Args:
        estadisticas_filtradas (list): Lista de diccionarios con estadísticas filtradas.
        ruta (str): Ruta del archivo de salida TSV.

    Returns:
        None
    """
    archivo = open(ruta, "w")

    # Escribir encabezado
    archivo.write("encabezado\tlongitud\tcontenido_gc\n")

    # Escribir cada secuencia que pasó los filtros
    for stats in estadisticas_filtradas:
        encabezado = stats["encabezado"]
        longitud = stats["longitud"]
        gc = stats["contenido_gc"]
        archivo.write(f"{encabezado}\t{longitud}\t{gc:.4f}\n")

    archivo.close()


def parsear_argumentos():
    """
    Parsea los argumentos de línea de comandos.

    Returns:
        argparse.Namespace: Objeto con los argumentos parseados.
    """
    parser = argparse.ArgumentParser(description="Analizador de secuencias FASTA")

    parser.add_argument("-i", "--input", required=True, help="Archivo FASTA de entrada")

    parser.add_argument("-o", "--output", required=True, help="Archivo TSV de salida")

    parser.add_argument("--min-len", type=int, default=0, help="Longitud mínima")

    parser.add_argument(
        "--max-len", type=int, default=float("inf"), help="Longitud máxima"
    )

    parser.add_argument("--min-gc", type=float, default=0, help="Contenido GC mínimo")

    parser.add_argument("--max-gc", type=float, default=1, help="Contenido GC máximo")

    args = parser.parse_args()

    return args


def main():
    """
    Función principal que coordina el procesamiento de archivos FASTA.

    Lee un archivo FASTA, calcula estadísticas, aplica filtros
    y guarda los resultados en un archivo TSV.
    """
    try:
        args = parsear_argumentos()

        print(f"Leyendo archivo: {args.input}")

        secuencias = leer_fasta(args.input)

        if not secuencias:
            print("Error: El archivo FASTA no contiene secuencias")
            sys.exit(1)

        print(f"  {len(secuencias)} secuencias encontradas")

        estadisticas = calcular_estadisticas(secuencias)

        estadisticas_filtradas = []

        for stats in estadisticas:
            cumple = pasa_filtros(
                stats, args.min_len, args.max_len, args.min_gc, args.max_gc
            )
            if cumple:
                estadisticas_filtradas.append(stats)

        print(f"  {len(estadisticas_filtradas)} secuencias pasan los filtros")

        escribir_resultados(estadisticas_filtradas, args.output)

        print(f"Resultados escritos en '{args.output}'")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()

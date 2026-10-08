"""
Создаём файл для работы с биологическими последовательностями.
Для этого реализуем два класса:
- Seq: хранение и анализ последовательности
- FastaReader: чтение файла в формате fasta
"""
class Seq:
    """
    Класс для хранения и анализа биологической последовательности.
    Хранит последовательность и её заголовок. Умеет определять алфавит (нуклеотидный или белковый) и считать GC-состав.

    Attributes
    ----------
    NUCLEOTIDES : set
        Множество символов нуклеотидов — A, C, G, T, U, N
    PROTEIN : set
        Множество символов 20 аминокислот
    _sequence : str
        Последовательность (защищённая)
    _header : str
        Заголовок FASTA-записи (защищённая)
    """
    NUCLEOTIDES=set("ACGTUN")
    PROTEIN=set("ACDEFGHIKLMNPQRSTVWY")

    def __init__(self, sequence, header=""):
        """
        Инициализирует объект последовательности.

        Parameters
        ----------
        sequence : str
            Последовательность (не пустая). Может быть в любом регистре —
            приведётся к верхнему.
        header : str, optional
            Заголовок FASTA-записи (по умолчанию пустая строка).

        Raises
        ------
        ValueError
            Если последовательность пустая
        """
        if not sequence:
            raise ValueError("Последовательность не может быть пустой")
        self._sequence=sequence.upper().strip()
        self._header=header.strip()

    @property
    def sequence(self):
        """
        Возвращает последовательность в верхнем регистре.

        Returns
        -------
        str
            Последовательность
        """
        return self._sequence
    @property
    def header(self):
        """
        Возвращает заголовок FASTA-записи.

        Returns
        -------
        str
            Заголовок
        """
        return self._header
    def __len__(self):
        """
        Возвращает длину последовательности.

        Returns
        -------
        int
            Количество символов
        """
        return len(self._sequence)
    def __str__(self):
        """
        Возвращает строковое представление объекта.

        Returns
        -------
        str
            Строка вида "Seq(header='...', length=N, alphabet=...)"
        """
        return f"Seq(header='{self._header}', length={len(self)}, alphabet={self.get_alphabet()})"
    def get_alphabet(self):
        """
        Определяет алфавит последовательности.

        Returns
        -------
        str
            Одно из значений:
                - "nucleotide" — если все символы ∈ NUCLEOTIDES
                - "protein" — если все символы ∈ PROTEIN
                - "unknown" — иначе
        """
        symbols= set(self._sequence)
        if symbols.issubset(self.NUCLEOTIDES):
            return "nucleotide"
        if symbols.issubset(self.PROTEIN):
            return "protein"
        return "unknown"
    def gc_content(self):
        """
        Считает GC-состав последовательности (доля G и C в %).

        Returns
        -------
        float
            GC-состав в процентах (0–100)

        Raises
        ------
        ValueError
            Если последовательность не нуклеотидная
        """
        if self.get_alphabet() !="nucleotide":
            raise ValueError("Только для нуклеотидов")
        g=self._sequence.count("G")
        c=self._sequence.count("C")
        return (g+c)/len(self)*100
class FastaReader:
    """
    Класс для чтения FASTA-файлов по записям.

    FASTA — текстовый формат для хранения последовательностей:
        >идентификатор описание
        ATGCGTACGT

    Читает файл как генератор — по одной записи за раз.
    Это позволяет обрабатывать файлы любого размера без загрузки
    всего содержимого в память.

    Attributes
    ----------
    _filename : str
        Путь к FASTA-файлу (защищённый) 
    """
    def __init__(self, filename):
        """
        Инициализирует читателя FASTA-файлов.

        Parameters
        ----------
        filename : str
            Путь к FASTA-файлу
        """
        self._filename=filename

    @property
    def filename(self):
        """
        Возвращает путь к файлу.

        Returns
        -------
        str
            Путь к FASTA-файлу
        """
        return self._filename
    def is_fasta(self):
        """
        Проверяет, что файл в формате FASTA.

        FASTA-файл должен начинаться с символа '>'.

        Returns
        -------
        bool
            True, если файл в формате FASTA, иначе False
        """
        try:
            with open(self._filename, encoding="utf-8") as file:
                first = file.readline().strip()
                return first.startswith(">")
        except FileNotFoundError:
            return False
    def read(self):
        """
        Генератор: читает файл по записям.

        Читает файл построчно, формирует записи и возвращает
        объекты Seq по одному. Не загружает весь файл в память.

        Yields
        ------
        Seq
            Объект последовательности для каждой FASTA-записи
        """
        header = None
        lines=[]
        with open(self._filename, encoding="utf-8") as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        seq="".join(lines)
                        if seq:
                            yield Seq(seq, header)
                    header = line[1:].strip()
                    lines=[]
                else:
                    lines.append(line)
        if header is not None:
            seq="".join(lines)
            if seq:
                yield Seq(seq,header)

# ДЕМОНСТРАЦИЯ
if __name__ == "__main__":
    # 1. Класс Seq
    seq = Seq("ATGCGTACGT", "ДНК")
    print(seq)
    print(f"Длина: {len(seq)}")
    print(f"Алфавит: {seq.get_alphabet()}")
    print(f"GC: {seq.gc_content():.2f}%")

    # 2. FastaReader
    reader = FastaReader("example.fasta")
    print(f"\nFASTA? {reader.is_fasta()}")

    for seq in reader.read():
        print(seq)
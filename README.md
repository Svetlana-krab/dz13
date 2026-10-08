# FASTA Reader

Минимальная библиотека для работы с биологическими последовательностями.

## Возможности

- Класс `Seq` — хранение и анализ последовательности:
  - определение алфавита (нуклеотидный/белковый)
  - расчёт GC-состава
  - длина последовательности
- Класс `FastaReader` — чтение FASTA-файлов:
  - проверка формата
  - чтение по записям (генератор — экономия памяти)

## Установка

Скопируйте файл `fasta_reader.py` в свой проект.

Требуется **Python 3.8+**. Внешних зависимостей нет.

## Использование

```python
from fasta_reader import Seq, FastaReader

# Класс Seq
seq = Seq("ATGCGTACGT", "ДНК")
print(len(seq))            # 10
print(seq.get_alphabet())  # nucleotide
print(seq.gc_content())    # 50.0

# Класс FastaReader
reader = FastaReader("example.fasta")
for s in reader.read():
    print(s)

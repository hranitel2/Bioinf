#!/bin/bash

# Создаём папку fastqs
mkdir -p fastqs

# Цикл от 1 до 10
for i in {1..10}; do
    # Имя файла
    filename="fastqs/sample_${i}.fastq"
    
    # Создаём файл и записываем строку с номером
    echo "This is sample number ${i}" > "$filename"
    
    echo "Created: $filename"
done

echo "All 10 samples generated in fastqs/"

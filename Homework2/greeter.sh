#!/bin/bash

# Проверяем, что передано ровно три аргумента
if [ $# -lt 3 ]; then
    echo "Error: Not enough arguments provided"
    exit 1
fi

# Запоминаем аргументы в читаемые переменные
FIRST_NAME="$1"
LAST_NAME="$2"
GROUP="$3"

# Выводим приветствие
echo "Welcome, $FIRST_NAME $LAST_NAME from group $GROUP!"

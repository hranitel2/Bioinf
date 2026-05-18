#!/bin/bash

# Проверяем, передан ли аргумент
if [ -z "$1" ]; then
    echo "No project name provided"
    exit 1
fi

# Запоминаем имя проекта
PROJECT_NAME="$1"

# Создаём папку проекта и подпапки
mkdir -p "$PROJECT_NAME"/{data,scripts,results}

# Создаём raw_data.txt и ставим права 600
touch "$PROJECT_NAME/data/raw_data.txt"
chmod 600 "$PROJECT_NAME/data/raw_data.txt"

# Создаём run_analysis.sh, делаем исполняемым и записываем скрипт
cat > "$PROJECT_NAME/scripts/run_analysis.sh" << 'EOF'
#!/bin/bash
echo "Hello from $1"
EOF

chmod +x "$PROJECT_NAME/scripts/run_analysis.sh"

# Финал
echo "Project '$PROJECT_NAME' initialized successfully."

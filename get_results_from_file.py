# Открываем файл для чтения
with open('bch(31,21)supercode(3_4).txt', 'r') as file:
    lines = file.readlines()

# Инициализируем списки для хранения значений
fer_list = []
ber_list = []

# Обрабатываем каждую строку
for line in lines:
    # Разделяем строку по запятым
    parts = line.strip().split(', ')
    
    # Извлекаем значения fer и ber
    fer = float(parts[0].split(' = ')[1])
    ber = float(parts[1].split(' = ')[1])
    
    # Добавляем значения в списки
    fer_list.append(fer)
    ber_list.append(ber)

# Выводим результаты
print("FER array:", fer_list)
print("BER array:", ber_list)
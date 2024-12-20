# Библиотека для обработки исходных данных сигналов акустических источников

**Библиотека для обработки исходных данных сигналов акустических источников** – набор инструментов, который позволяет проводить предварительную обработку и суммирование исходных акустических распределённых сигналов (в том числе от пневматических источников). Программные компоненты модуля позволяют различными методами провести предварительную фильтрацию исходных распределенных сигналов и рассчитать итоговый сигнал в заданной точке пространства, с учетом его затухания и отражения.

![Logotype](./docs/НТИ_логотип_с_плашкой_RGB.png)

## Установка

Установка не предполагается.

<!---
## Разработка

Разработка UI приветствуется.
-->

### Сборка

Используется система сборки CMake.

Главный файл находится в каталоге `app\gund_solvers`.

Для сборки тестов указывается опция `WITH_TESTS`.

Используемые сторонние библиотеки:
- [nlohmann](https://github.com/nlohmann/json/releases/)
- [Simple-FFT](https://github.com/d1vanov/Simple-FFT)


## Функциональность

Учтено отражение от водной поверхности.

Исходные акустические распределённые сигналы зачитываются из файлов `*.sig`.

Поддерживаются различные варианты предварительной фильтрации исходных распределенных сигналов.

## Благодарности

Работа выполнена Инжиниринговым центром по трудноизвлекаемым полезным ископаемым Центра компетенций НТИ на базе МФТИ по направлению "Искусственный интеллект" в рамках "дорожной карты" развития высокотехнологичного направления "Искусственный интеллект" на период до 2030 года при поддержке Фонда НТИ.

<!---
## Авторы
-->

## Contributing

Комментарии и улучшения приветствуются.

<!---
## Links

- Project homepage: 
- Repository: 
- Related projects:
-->

## Лицензия

Распространяется по лицензии BSD 3-Clause.

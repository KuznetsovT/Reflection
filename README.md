# Reflection

Автоматизация поиска позиций гониомера четырехкружного дифрактометра с каппа-геометрией для измерения дифракционного рефлекса кристалла с заданной матрицей ориентации.


Лабораторные монокристальные рентгеновские дифрактометры, оснащенные многокружными гониометрами и 2D-детекторами, являются основными современными приборами для проведения рентгеноструктурного анализа.
Распространенные программы для анализа дифракционных(Брэгговских) пиков заточены на определение матрицы ориентации и проверку нахождения заданного рефлекса в отражающем положении по заданной позиции. 
В данном докладе рассматривается обратная задача: нахождение позиций гониометра для заданного рефлекса, при которых данный рефлекс находится в отражающем положении. Предлагаемое решение исключает перебор всех возможных позиций гониометра с некоторым шагом. 
Данное решение можно разбить на этапы:
1. Нахождение решения для одноосного отражателя (кристалл может двигаться только вокруг одной оси, перпендикулярной направлению падающего излучения).
2. Нахождение частного решения для гониометра, подчиняющегося Эйлеровой геометрии, удовлетворяющего определённым заданным свойствам. 
3. Нахождение и систематизация множества решений для гониометра с Эйлеровой геометрией. 
4. Автоматизация поиска путём написания на С++ соответствующей программы. 
Планируется дальнейшее расширение функциональности лабораторного дифрактометра, уточнением встроенной модели столкновений, для возможности планирования съёмки, не предусмотренной стандартным ПО. Например — автоматический поиск стратегии измерения отдельных высокоугловых рефлексов (с расщеплением Ka1/Ka2) для прецизионного определения параметров элементарных ячеек.


Для знакомвства с диффракцинной физикой посмотрите файл Diffraction.h.
В этом файле реализуется условие Вульфа-Брэгга (можно поискать про это в интернете), описывающее условие диффракции в кристалле.
Также в файле содержится информация об индексах Миллера, что необходимо для нахождения векторов обратного пространства(с которыми в дальнейшем очень удобно работать)

Так как большинство вычислений векторные, в классе R3 содержатся базовые элементы для векторной алгебры.

Класс Gonio непосредственно отвечает за предсказание позиций гониометра. В нем есть метод "omega_rotation", который находит решение для одноосного отражателя.
Далее, чтобы систематизировать решения трехосного отражателя, введем две координаты - кси и пси.
Кси - это угол между рефлексом и горизонтальной осью. Примерно как геометрический круг.


Тогда найдётся ось, вращение вокруг которой не будет менять положение рефлекса, а значит угол поворота вокруг этой оси мы назовём Пси.
Пси=0 когда угол хи(углол поворота гониометра в Эйлеровой геометрии) минимален.

Теперь смотрим на функцию "Gonio::diff_rotation". Именно эта функция является трехосным предсказателем. 
Кроме того, результат может выводиться как в Эйлеровой геометрии, так и в нативной каппа-геометрии.
Структура Solution нам дополнительно говорит, сколько у нас получается решений, потому что если хи равен нулю, параметры омега и фи линейно зависимы.

Дополнительно, недавно появился простенький предсказатель столкновений, который опирается на список проверенных с некоторым шагом позиций.
     
     
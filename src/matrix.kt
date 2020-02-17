import java.io.File
import java.io.InputStream

data class GaussParam(val m: Matrix, val e: Matrix, val D: Double, val i: Int, val j: Int)


class Matrix (val height: Int, val width: Int) {

    private val data = Array(height) { DoubleArray(width) { 0.0 } }

    fun get(row: Int, column: Int): Double {
        if (row in 0 until height && row in 0 until width) {
            return data[row][column]
        }
        return 0.0
    }

    fun set(row: Int, column: Int, value: Double) {
        if (row in 0 until height && row in 0 until width) {
            data[row][column] = value
        }
    }

    fun print(): String {
        var to = ""
        for (row in data) {
            for (x in row) {
                print("$x ")
                to = to + "$x "
            }
            println()
            to = to + "\n"
        }
        return to
    }

    fun mul(scal: Double): Matrix {
        var m = Matrix(height, width)
        for (i in 0 until height) {
            for (j in 0 until width) {
                m.set(i, j,data[i][j] * scal)
            }
        }
        return m
    }

    fun add(rhs: Matrix): Matrix? {
        if (height == rhs.height && width == rhs.width) {
            var m = Matrix(height, width)
            for (i in 0 until height) {
                for (j in 0 until width) {
                    m.set(i, j, data[i][j] + rhs.get(i, j))
                }
            }
            return m
        }
        return null
    }

    fun transp(): Matrix {
        var m = Matrix(width, height)
        for (i in 0 until m.height) {
            for (j in 0 until m.width) {
                m.set(i, j, data[j][i])
            }
        }
        return m
    }

    fun mul(rhs: Matrix): Matrix? {
        if (width == rhs.height) {
            var m = Matrix(width, height)
            for (i in 0 until m.height) {
                for (j in 0 until m.width) {
                    for (k in 0 until width)
                        m.set(i, j, m.get(i, j) + data[i][k] * rhs.get(k, j))
                }
            }
            return m
        }
        return null
    }

    private fun step1(vals: GaussParam): GaussParam {
        var m = vals.m
        var e = vals.e
        var D = vals.D
        var i = vals.i
        var j = vals.j
        var tmp = 0.0
        var inv = false
        if (height == width) {
            inv = true
        }
        while(true) {
            if (m.get(i, j) == 0.0) {
                step2(GaussParam(m, e, D, i, j))
                return GaussParam(m, e, D, i, j)
            }
            if (i == j) {
                D = D * m.get(i, j)
            }
            tmp = m.get(i, j)
            for (t in 0 until width) {
                m.set(i, t, m.get(i, t) / tmp)
                if (inv) {
                    e.set(i, t, e.get(i, t) / tmp)
                }
            }
            if (i == height - 1) {
                return GaussParam(m, e, D, i, j)
            }
            for (t in i + 1 until height) {
                tmp = m.get(t, j)
                for (k in 0 until width) {
                    m.set(t, k, m.get(t, k) - m.get(i, k) * tmp)
                    if (inv) {
                        e.set(t, k, m.get(t, k) - m.get(i, k) * tmp)
                    }
                }
            }
            if (j == width - 1) {
                return GaussParam(m, e, D, i, j)
            }
            i++
            j++
        }
    }

    private fun step2(vals: GaussParam) {
        var m = vals.m
        var e = vals.e
        var D = vals.D
        var i = vals.i
        var j = vals.j
        var tmp = 0.0
        var inv = false
        if (height == width) {
            inv = true
        }

        if (i < height - 1) {
            for (t in i + 1 until height - 1) {
                if (m.get(t, j) != 0.0) {
                    for (k in 0 until width) {
                        var temp = m.get(i, k)
                        m.set(i, k, m.get(t, k))
                        m.set(t, k, temp)
                        if (inv) {
                            var temp = e.get(i, k)
                            e.set(i, k, e.get(t, k))
                            e.set(t, k, temp)
                        }
                    }
                    D = D * -1;
                    step1(GaussParam(m, e, D, i, j))
                    return
                }
            }
        }
        D = 0.0;
        if (j == width - 1) {
            return
        }
        ++j;
        step1(GaussParam(m, e, D, i, j))
    }

    private fun gausselim(): Triple<Matrix, Matrix, Double> {
        var m = Matrix(height, width)
        for (i in 0 until height) {
            for (j in 0 until width) {
                m.set(i, j, data[i][j])
            }
        }
        var e = Matrix(height, width)
        for (i in 0 until height) {
            e.set(i, i, 1.0)

        }
        var D: Double = 1.0
        var i = 0
        var j = 0
        var vals = step1(GaussParam(m, e, D, i, j))
        return Triple(vals.m, vals.e, vals.D)
    }

    fun det(): Double? {
        if (height == width) {
            var vals = gausselim()
            return vals.third
        }
        return null
    }

    fun inverse(): Matrix? {
        if (height == width) {
            var vals = gausselim()
            return vals.second
        }
        return null
    }

    fun Gauss(): Matrix {
        var vals = gausselim()
        return vals.first
    }

    fun rank(): Int {
        var vals = gausselim()
        return vals.first.height
    }
}

fun main() {
    var matrices = arrayListOf<Matrix>()

    val inputStream: InputStream = File("matrix.txt").inputStream()
    val input = mutableListOf<String>()
    inputStream.bufferedReader().useLines { lines -> lines.forEach { input.add(it)} }

    var it = 0
    for(k in it until input.size) {
        var last = 0
        if (k == input.size - 1) {
            last = 1
        }
        if (input[k] == "" || last == 1) {
            var line = input[it]
            var nums = line.split("\t")
            var m = Matrix(k + last - it, nums.size)
            for (i in 0 until m.height) {
                line = input[it++]
                nums = line.split(" ")
                for (j in 0 until m.width) {
                    if (nums[j] != "") {
                        m.set(i, j, nums[j].toDouble())
                    }
                }
            }
            matrices.add(m)
            it = k + 1
        }
    }

    var ops = mutableMapOf<Char, String>()
    ops['a'] = "addition"
    ops['b'] = "matrix multiplication"
    ops['c'] = "scalar multiplication"
    ops['d'] = "transposition"
    ops['e'] = "inverse"
    ops['f'] = "rank"
    ops['g'] = "determinant"
    ops['h'] = "Gauss elimination"
    ops['i'] = "Other: new matrix"

    var result = Matrix(0,0)

    while(true) {
        println("Current matrices:")
        matrices.forEachIndexed { i, m ->
            println("${i + 1}.")
            m.print()
            println()
        }
        println()
        println("Possible operations:")
        for(op in ops) {
            var key = op.key
            var value = op.value
            if(key == 'i') {
                println("Other:")
            }
            println("\t$key) $value")
        }

        println()
        println("Choose an operation!")

        var c = readLine()
        c?.toLowerCase()
        val op = c?.get(0)
        if (op in 'a'..'i') {
            if (op in 'a'..'e') {
                println("Save result to the old/new matrix [o/n]")
                c = readLine()
                c?.toLowerCase()
            }
            val save = c?.get(0)
            var m = 0
            println()
            matrices.forEachIndexed { i, m ->
                println("${i + 1}.")
                m.print()
                println()
            }
            if (op in 'a'..'b') {
                println("Choose two of the matrices above!")
                c = readLine()
                m = c!!.toInt() - 1
                c = readLine()
                val m2 = c!!.toInt() - 1
                if(op == 'a') {
                    result = matrices[m].add(matrices[m2])!!
                }
                if(op == 'b') {
                    result = matrices[m].mul(matrices[m2])!!
                }
            }
            else if (op != 'i') {
                    println("Choose one of the matrices above!")
                    c = readLine()
                    m = c!!.toInt() - 1
                when (op) {
                    'c' -> {
                        println("Enter a scalar!")
                        c = readLine()
                        val scal = c!!.toDouble()
                        result = matrices[m].mul(scal)
                    }
                    'd' -> {
                        result = matrices[m].transp()
                    }
                    'e' -> {
                        result = matrices[m].inverse()!!
                    }
                    'f' -> {
                        println(matrices[m].rank())
                    }
                    'g' -> {
                        println(matrices[m].det())
                    }
                    'h' -> {
                        result = matrices[m].Gauss()
                    }
                }
            }
            else {
                println("Enter height and width!")
                c = readLine()
                val h = c!!.toInt()
                c = readLine()
                val w = c!!.toInt()

                var new = Matrix(h, w)
                for (i in 0 until new.height) {
                    var line = readLine()
                    var nums = line!!.split("\t")
                    for (j in 0 until new.width) {
                        if (nums[j] != "") {
                            new.set(i, j, nums[j].toDouble())
                        }
                    }
                }
                matrices.add(new)
            }

            var new = result

            if(save == 'o') {
                matrices[m] = result
            }
            if(save == 'n') {
                matrices.add(result)
            }

            var to = ""

            println()
            println()
            for (m in matrices) {
                to = to + m.print() + "\n"
                println()
            }

            File("matrix.txt").bufferedWriter().use { out -> out.write(to) }

            println("Do you want to exit? [y/n]")
            c = readLine()
            c?.toLowerCase()
            var e = c?.get(0)
            if(e == 'y') {
                return
            }


        }

    }


}




-- Copyright (c) 2025 James Cook
-- Matrix calculation functions

-- Matrix functions:
--
-- global type matrix(sequence s)
-- global function NewMatrix(integer rows, integer cols)
-- global function GetMatrixRows(sequence a)
-- global function GetMatrixCols(sequence a)
-- global function MatrixMultiply(matrix a, matrix b) -- MatrixMultiplication
-- global function MatrixTransformation(matrix a)


global function GetMatrixRows(sequence a)
    return length(a)
end function

global function GetMatrixCols(sequence a)
    if length(a) then
        return length(a[1])
    end if
    return 0
end function

global type matrix(sequence s)
    integer lenRows, lenCols
    lenRows = GetMatrixRows(s)
    if lenRows = 0 then
        return 0
    end if
    lenCols = GetMatrixCols(s)
    for i = 2 to lenRows do
        if length(s[i]) != lenCols then
            return 0
        end if
    end for
    return 1
end type

global function NewMatrix(integer rows, integer cols)
    return repeat(repeat(0, cols), rows)
end function

global function MatrixMultiply(matrix a, matrix b)
-- ret[i] =
-- {
--  a[i][k1] * b[k1][j1] + a[i][k2] * b[k2][j1],
--  a[i][k1] * b[k1][j2] + a[i][k2] * b[k2][j2],
--  a[i][k1] * b[k1][j3] + a[i][k2] * b[k2][j3],
--  a[i][k1] * b[k1][j4] + a[i][k2] * b[k2][j4]
-- }
-- row0 = ret[i]
-- row1 = a[i]
    integer rows, cols, len
    sequence row0, row1, sum, ret
    matrix ret
    len = GetMatrixRows(b)
    if GetMatrixCols(a) != len then
        abort(1)
    end if
    rows = GetMatrixRows(a)
    cols = GetMatrixCols(b)
    ret = NewMatrix(rows, cols)
    for i = 1 to rows do -- rows of "a"
        row0 = ret[i]
        row1 = a[i]
        for j = 1 to cols do -- cols of "b"
            sum = 0
            for k = 1 to len do -- k is cols of "a", rows of "b"
                sum += row1[k] * b[k][j]
            end for
            row0[j] = sum
        end for
        ret[i] = row0
    end for
    return ret
end function

global function MatrixTransformation(matrix a)
    integer rows, cols
    sequence ret, tmp
    rows = GetMatrixRows(a)
    cols = GetMatrixCols(a)
    ret = NewMatrix(cols, rows)
    for row = 1 to rows do
        tmp = a[row]
        for col = 1 to cols do
            ret[col][row] = tmp[col]
        end for
    end for
    return ret
end function

-- See also:
-- https://www.purplemath.com/modules/mtrxmult3.htm

-- end of file.

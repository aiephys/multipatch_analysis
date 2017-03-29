# *-* coding: utf-8 *-*
import pyqtgraph as pg


class MatrixItem(pg.QtGui.QGraphicsItemGroup):
    """GraphicsItem displaying a table with column / row labels and text in
    each cell.
    
    Parameters
    ----------
    text : 2d array or nested lists
        Strings to display inside each cell
    fgcolor : 2d array or nested lists
        Text colors for each cell
    bgcolor : 2d array or nested lists
        Background colors for each cell
    rows : 1d array or list
        Strings to display as row header
    cols : 1d array or list
        Strings to display as col header
    size : float
        Width of each cell
    """
    def __init__(self, text, fgcolor, bgcolor, rows, cols, size=50):
        pg.QtGui.QGraphicsItemGroup.__init__(self)

        for i,row in enumerate(rows):
            txt = pg.QtGui.QGraphicsTextItem(row, parent=self)
            br = txt.boundingRect()
            txt.setPos(-br.width() - 10, i * size + size/2. - br.center().y())
            txt.setDefaultTextColor(pg.mkColor('w'))

        for i,col in enumerate(cols):
            txt = pg.QtGui.QGraphicsTextItem(col, parent=self)
            br = txt.boundingRect()
            txt.setPos(i * size + size/2 - br.center().x(), -br.height() - 10)
            txt.setDefaultTextColor(pg.mkColor('w'))

        for i,row in enumerate(rows):
            for j,col in enumerate(cols):
                x = j * size
                y = i * size
                rect = pg.QtGui.QGraphicsRectItem(x, y, size, size, parent=self)
                rect.setBrush(pg.mkBrush(bgcolor[i][j]))
                rect.setZValue(-10)

                txt = pg.QtGui.QGraphicsTextItem(text[i][j], parent=self)
                br = txt.boundingRect()
                txt.setPos(x + size/2 - br.center().x(), y + size/2 - br.center().y())
                txt.setDefaultTextColor(pg.mkColor(fgcolor[i][j]))

        br = pg.QtCore.QRectF()
        for item in self.childItems():
            br = br.united(self.mapRectFromItem(item, item.boundingRect()))
        self._bounding_rect = br

    def boundingRect(self):
        return self._bounding_rect
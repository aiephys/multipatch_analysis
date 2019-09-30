import time
from .interactive_mode import interactive_mode


class ProgressBar:
    """A progress bar that works nicely in any user interface (or lack thereof).
    
    Supports Qt, Jupyter notebook, console, and headless modes.
    
    Example::
    
        with ProgressBar(message, maximum) as prg_bar:
            # do some work here
            try:
                prg_bar.update(new_value, status_message)
            except prg_bar.CanceledError:
                # handle user cancellation here
                
    """
    class CanceledError(Exception):
        pass
    
    def __init__(self, message, maximum, update_interval=0.1, mode=None):
        self.message = message
        self.maximum = maximum
        self.value = 0

        if mode is None:
            mode = interactive_mode()
        
        self.mode = mode
        self.update_dt = update_interval
        self._last_update = 0
        self._need_update = True

        if self.mode == 'qt':
            import pyqtgraph as pg
            self.dlg = pg.ProgressDialog(self.message, maximum=1000)
        elif self.mode == 'file':
            print(message)
            # logging to file; don't need frequent updates.
            self.update_dt = 30.0
        else:
            print(message)
            self.last_line_len = 0

    def update(self, value, status):
        now = time.time()
        self.value = value
        self._need_update = True

        # rate-limit updates
        if value < self.maximum and now < self._last_update + self.update_dt:
            return
        self._last_update = now
        self._need_update = False
        
        if self.mode == 'qt':
            self.dlg.setLabelText(self.message + '\n' + status)
            self.dlg.setValue(1000. * value / self.maximum)
            if self.dlg.wasCanceled():
                raise ProgressBar.CanceledError()
        elif self.mode in ('tty', 'ipynb'):
            nh = int(20.* value / self.maximum)
            complete = "#" * nh
            incomplete = " " * (20 - nh)
            if self.mode == 'ipynb':
                import IPython.display
                IPython.display.clear_output(wait=True)
                print(self.message)
            line = '  [%s%s]  %s' % (complete, incomplete, status)
            line_len = len(line.split('\n')[-1])
            spaces = ' ' * max(0, self.last_line_len - line_len)
            self.last_line_len = line_len
            print('\r' + line + spaces, end='')
        else:
            print('  ' + status)

        if value >= self.maximum and self.mode != 'qt':
            print("\n  done.")

    def __enter__(self):
        return self
        
    def __exit__(self, exType, exValue, exTrace):
        if self.mode == 'qt':
            self.dlg.setValue(self.dlg.maximum())
            if self.dlg.isVisible():
                self.dlg.close()

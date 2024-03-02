'''zstd interface via the open() command'''
## Thanks to John Do and jasonharper
## https://stackoverflow.com/questions/64488350/\
## using-zstandard-to-compress-a-file-in-python

import builtins
import io
import zstandard as zstd

class ZstdReader:
    '''zstd reader class for "with open(...) as ..." context manager '''
    def __init__(self, filename,mode):
        self.filename = filename
        if 't' in mode and 'b' in mode:
            raise ValueError(f'Invalid mode={mode}; has both "t" and "b"')
        self.binarymode = bool('b' in mode)
        self.fp = None
        self.bin_reader = None
        self.txt_reader = None

    def __enter__(self):
        self.fp = builtins.open(self.filename, 'rb')
        dctx = zstd.ZstdDecompressor()
        self.bin_reader = dctx.stream_reader(self.fp)
        if self.binarymode:
            return self.bin_reader
        self.txt_reader = io.TextIOWrapper(self.bin_reader,
                                           encoding='utf-8')
        return self.txt_reader

    def __exit__(self, *args):
        self.fp.close()
        return False

class ZstdWriter:
    '''zstd writer class for "with open(...) as ..." context manager '''
    def __init__(self, filename,mode):
        self.filename = filename
        if 't' in mode and 'b' in mode:
            raise ValueError(f'Invalid mode={mode}; has both "t" and "b"')
        self.binarymode = bool('b' in mode)
        self.fp = None
        self.bin_writer = None
        self.txt_writer = None

    def __enter__(self):
        self.fp = builtins.open(self.filename, 'wb')
        ctx = zstd.ZstdCompressor()
        self.bin_writer = ctx.stream_writer(self.fp)
        if self.binarymode:
            return self.bin_writer
        self.txt_writer = io.TextIOWrapper(self.bin_writer,
                                           encoding='utf-8')
        return self.txt_writer

    def __exit__(self, *args):
        if self.binarymode:
            self.txt_writer.flush()
        self.bin_writer.flush(zstd.FLUSH_FRAME)
        self.fp.close()
        return False

def open(filename, mode='r'): #pylint: disable=redefined-builtin
    '''like builtin.open but for zst files'''
    if 'w' in mode:
        return ZstdWriter(filename,mode)
    return ZstdReader(filename,mode)

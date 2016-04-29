from stammp.utils import execute

class PipelineModule:

    def __init__(self):
        self._data = {}

    def prepare(self):
        pass

    def execute(self):
        pass

    def cleanup(self):
        pass

    def get(self, key='output', default=None):
        return self._data.get(key, default)


class CmdPipelineModule(PipelineModule):

    def __init__(self, remove_files=True):
        super().__init__()
        self._remove_files = remove_files
        self._data['files'] = []
        self._cmds = []

    def execute(self):
        for cmd in self._cmds:
            yield execute(cmd)

    def cleanup(self):
        if self._remove_files:
            for rm_file in self._data['files']:
                execute('rm -rf %s' % rm_file, exit=False)


class DummyModule(PipelineModule):

    def __init__(self, infile):
        super().__init__()
        self._data['output'] = infile


class Pipeline:

    def __init__(self, initial_file):
        dmod = DummyModule(initial_file)
        self._jobs = [dmod]
        self._current = 0

    def schedule(self, module):
        self._jobs.append(module)

    def __iter__(self):
        self._current = 0
        return self

    def __next__(self):
        if self._current < len(self._jobs):
            self._current += 1
            return self._jobs[self._current - 1]
        else:
            raise StopIteration

    def cur_output(self):
        for job in self._jobs[::-1]:
            if job.get() is not None:
                return job.get()
        return None

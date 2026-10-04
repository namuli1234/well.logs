import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backends.backend_pdf import PdfPages
import petro as P

DISCLAIMER = ('Interpretations are engineering estimates that depend on data quality, environmental '
              'corrections, model assumptions and parameter choice. Review all results.')
PARAMS = [('GR min', 'grmin', '15'), ('GR max', 'grmax', '120'), ('Matrix density (g/cc)', 'rma', '2.65'),
          ('Fluid density (g/cc)', 'rfl', '1.0'), ('Rw (ohm-m)', 'rw', '0.08'), ('Rsh (ohm-m)', 'rsh', '2.0'),
          ('Archie a', 'a', '1.0'), ('Archie m', 'm', '2.0'), ('Archie n', 'n', '2.0'),
          ('Vsh cutoff (<)', 'cvsh', '0.35'), ('Porosity cutoff (>)', 'cphi', '0.10'), ('Sw cutoff (<)', 'csw', '0.60')]


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('WellLog Pro - Petrophysical Interpretation')
        self.geometry('1400x850')
        self.w = self.mapped = self.res = None
        side = ttk.Frame(self, padding=6); side.pack(side='left', fill='y')
        ttk.Button(side, text='Load LAS / CSV...', command=self.load).pack(fill='x', pady=2)
        self.v = {}
        for lab, key, dflt in PARAMS:
            r = ttk.Frame(side); r.pack(fill='x')
            ttk.Label(r, text=lab, width=22).pack(side='left')
            self.v[key] = tk.StringVar(value=dflt)
            ttk.Entry(r, textvariable=self.v[key], width=8).pack(side='right')
        self.opt = {}
        for lab, key, vals in [('Vsh method', 'vsh', ['Linear', 'Larionov (tertiary)', 'Larionov (older)', 'Clavier', 'Steiber']),
                               ('Porosity', 'por', ['Density', 'Neutron', 'Density-Neutron average']),
                               ('Sw model', 'sw', ['Archie', 'Simandoux', 'Indonesia'])]:
            ttk.Label(side, text=lab).pack(anchor='w')
            self.opt[key] = ttk.Combobox(side, values=vals, state='readonly'); self.opt[key].current(0)
            self.opt[key].pack(fill='x')
        ttk.Label(side, text='Matrix preset').pack(anchor='w')
        mp = ttk.Combobox(side, values=list(P.MATRIX), state='readonly'); mp.pack(fill='x')
        mp.bind('<<ComboboxSelected>>', lambda e: self.v['rma'].set(str(P.MATRIX[mp.get()])))
        ttk.Label(side, text='Zones (Name,top,base per line)').pack(anchor='w', pady=(6, 0))
        self.zones = tk.Text(side, width=30, height=5); self.zones.pack(fill='x')
        ttk.Button(side, text='CALCULATE', command=self.calculate).pack(fill='x', pady=6)
        ttk.Button(side, text='Export CSV', command=self.exp_csv).pack(fill='x', pady=1)
        ttk.Button(side, text='Export PNG (current tab)', command=self.exp_png).pack(fill='x', pady=1)
        ttk.Button(side, text='Export PDF report', command=self.exp_pdf).pack(fill='x', pady=1)
        ttk.Label(side, text=DISCLAIMER, wraplength=230, foreground='#a00', font=('Segoe UI', 8)).pack(pady=8)

        self.nb = ttk.Notebook(self); self.nb.pack(side='right', fill='both', expand=True)
        self.figs, self.canv = {}, {}
        for tab in ('Log Viewer', 'Crossplots'):
            f = ttk.Frame(self.nb); self.nb.add(f, text=tab)
            fig = Figure(figsize=(10, 7), dpi=100); c = FigureCanvasTkAgg(fig, f)
            NavigationToolbar2Tk(c, f); c.get_tk_widget().pack(fill='both', expand=True)
            self.figs[tab], self.canv[tab] = fig, c
        self.txt = {}
        for tab in ('QC', 'Zones / Summary'):
            f = ttk.Frame(self.nb); self.nb.add(f, text=tab)
            t = tk.Text(f, font=('Consolas', 10)); t.pack(fill='both', expand=True); self.txt[tab] = t
        self.cross = []

    # ---------- helpers ----------
    def prm(self):
        return {k: float(s.get()) for k, s in self.v.items()}

    def cur(self, std):
        n = self.mapped.get(std)
        return self.w['data'][n] if n else None

    def settext(self, tab, s):
        self.txt[tab].delete('1.0', 'end'); self.txt[tab].insert('end', s)

    # ---------- load ----------
    def load(self):
        p = filedialog.askopenfilename(filetypes=[('Well logs', '*.las *.LAS *.csv'), ('All', '*.*')])
        if not p: return
        try:
            self.w = P.parse_csv(p) if p.lower().endswith('.csv') else P.parse_las(p)
            self.mapped = P.map_curves(self.w['names'])
        except Exception as e:
            messagebox.showerror('Load failed', str(e)); return
        self.res = None
        self.settext('QC', P.qc_report(self.w, self.mapped))
        gr = self.cur('GR')
        if gr is not None and np.isfinite(gr).any():   # auto GR baselines: 5th/95th percentile, editable
            self.v['grmin'].set(f'{np.nanpercentile(gr, 5):.1f}'); self.v['grmax'].set(f'{np.nanpercentile(gr, 95):.1f}')
        d = self.w['depth']
        self.zones.delete('1.0', 'end'); self.zones.insert('end', f'Whole interval,{np.nanmin(d):.1f},{np.nanmax(d):.1f}')
        self.plot_logs(); self.nb.select(0)

    # ---------- calculate ----------
    def calculate(self):
        if not self.w: messagebox.showinfo('No data', 'Load a LAS/CSV file first.'); return
        try: p = self.prm()
        except ValueError: messagebox.showerror('Parameters', 'All parameters must be numbers.'); return
        gr, rt, rhob, nphi = self.cur('GR'), self.cur('RT'), self.cur('RHOB'), self.cur('NPHI')
        miss = [n for n, c in (('GR', gr), ('Deep resistivity (RT)', rt)) if c is None]
        pm = self.opt['por'].get()
        if pm != 'Neutron' and rhob is None: miss.append('RHOB')
        if pm != 'Density' and nphi is None: miss.append('NPHI')
        if miss: messagebox.showerror('Missing curves', 'Required curves not found: ' + ', '.join(miss) + '\nCurves are never assumed.'); return
        v = P.vsh(gr, p['grmin'], p['grmax'], self.opt['vsh'].get())
        pd_ = P.phi_density(rhob, p['rma'], p['rfl']) if rhob is not None else None
        pn = np.clip(P.norm_nphi(nphi), 0, 0.6) if nphi is not None else None
        phi = pd_ if pm == 'Density' else pn if pm == 'Neutron' else (pd_ + pn) / 2
        sm = self.opt['sw'].get()
        if sm == 'Archie': sw = P.sw_archie(rt, phi, p['rw'], p['a'], p['m'], p['n'])
        elif sm == 'Simandoux': sw = P.sw_simandoux(rt, phi, v, p['rw'], p['rsh'], p['a'], p['m'])
        else: sw = P.sw_indonesia(rt, phi, v, p['rw'], p['rsh'], p['a'], p['m'], p['n'])
        pay = P.pay_flag(v, phi, sw, {'vsh': p['cvsh'], 'phi': p['cphi'], 'sw': p['csw']})
        self.res = dict(VSH=v, PHI=phi, SW=sw, SH=1 - sw, PAY=pay)
        self.used = dict(p, vsh_method=self.opt['vsh'].get(), por_method=pm, sw_model=sm)
        self.plot_logs(); self.plot_cross(); self.summary(); self.nb.select(0)

    # ---------- plots ----------
    def plot_logs(self):
        fig = self.figs['Log Viewer']; fig.clf(); d = self.w['depth']
        tracks = []
        if self.cur('GR') is not None: tracks.append(('GR', [('GR', 'green', None)]))
        res = [(s, c) for s, c in (('RT', 'red'), ('RS', 'blue'), ('RXO', 'black')) if self.cur(s) is not None]
        if res: tracks.append(('Resistivity', [(s, c, None) for s, c in res]))
        if self.cur('RHOB') is not None or self.cur('NPHI') is not None: tracks.append(('RHOB / NPHI', None))
        if self.cur('DT') is not None: tracks.append(('DT', [('DT', 'purple', None)]))
        if self.res: tracks.append(('Vsh / PHI / Sw (calculated)', None))
        if not tracks: fig.text(.3, .5, 'No recognised curves to plot'); self.canv['Log Viewer'].draw(); return
        axs = fig.subplots(1, len(tracks), sharey=True); axs = np.atleast_1d(axs); self.vl = []
        for ax, (title, specs) in zip(axs, tracks):
            ax.set_title(title, fontsize=8); ax.grid(alpha=.3)
            if title == 'RHOB / NPHI':
                if self.cur('RHOB') is not None:
                    ax.plot(self.cur('RHOB'), d, 'r', lw=.8); ax.set_xlim(1.95, 2.95); ax.set_xlabel('RHOB', color='r', fontsize=7)
                if self.cur('NPHI') is not None:
                    a2 = ax.twiny(); a2.plot(P.norm_nphi(self.cur('NPHI')), d, 'b', lw=.8); a2.set_xlim(.45, -.15)
                    a2.set_xlabel('NPHI', color='b', fontsize=7)
            elif title.startswith('Vsh'):
                for k, c in (('VSH', 'brown'), ('PHI', 'blue'), ('SW', 'cyan')): ax.plot(self.res[k], d, c, lw=.9, label=k)
                ax.fill_betweenx(d, 0, 1, where=self.res['PAY'], color='gold', alpha=.35, label='PAY')
                ax.set_xlim(0, 1); ax.legend(fontsize=6, loc='lower right')
            else:
                for s, c, _ in specs: ax.plot(self.cur(s), d, c, lw=.8, label=s)
                if title == 'Resistivity': ax.set_xscale('log'); ax.set_xlim(.2, 2000); ax.legend(fontsize=6)
            self.vl.append(ax.axhline(d[0], color='k', lw=.6, visible=False))
        axs[0].set_ylim(np.nanmax(d), np.nanmin(d)); axs[0].set_ylabel(f"Depth ({self.w['depth_unit']})")
        fig.suptitle(f"{self.w['well']}   [measured curves vs. calculated (right track)]", fontsize=9)
        fig.tight_layout(); c = self.canv['Log Viewer']
        def mv(e):   # synchronised depth cursor
            if e.inaxes is None or e.ydata is None: return
            for l in self.vl: l.set_ydata([e.ydata] * 2); l.set_visible(True)
            c.draw_idle()
        fig.canvas.mpl_connect('motion_notify_event', mv); c.draw()

    def plot_cross(self):
        fig = self.figs['Crossplots']; fig.clf(); r = self.res; p = self.used
        rhob, nphi, rt, d = self.cur('RHOB'), self.cur('NPHI'), self.cur('RT'), self.w['depth']
        ax1, ax2 = fig.subplots(1, 2)
        if rhob is not None and nphi is not None:
            s = ax1.scatter(P.norm_nphi(nphi), rhob, c=self.cur('GR'), s=4, cmap='viridis')
            fig.colorbar(s, ax=ax1, label='GR'); ax1.set_xlim(-.05, .5); ax1.set_ylim(3.0, 1.9)
            ax1.set_xlabel('NPHI (v/v)'); ax1.set_ylabel('RHOB (g/cc)')
        ax1.set_title('Density-Neutron crossplot', fontsize=9); ax1.grid(alpha=.3)
        ok = np.isfinite(r['PHI']) & np.isfinite(rt) & (r['PHI'] > 0) & (rt > 0)
        s = ax2.scatter(r['PHI'][ok], rt[ok], c=d[ok], s=4, cmap='plasma'); fig.colorbar(s, ax=ax2, label='Depth')
        x = np.logspace(-2, 0, 50)
        for sw in (1.0, .75, .5, .25):
            ax2.plot(x, p['a'] * p['rw'] / (x ** p['m'] * sw ** p['n']), 'k--', lw=.7); 
            ax2.text(.5, p['a'] * p['rw'] / (.5 ** p['m'] * sw ** p['n']), f'Sw={sw:g}', fontsize=7)
        ax2.set_xscale('log'); ax2.set_yscale('log'); ax2.set_xlim(.01, 1); ax2.set_ylim(.1, 1000)
        ax2.set_xlabel('PHI (v/v)'); ax2.set_ylabel('Rt (ohm-m)')
        ax2.set_title(f"Pickett plot (a={p['a']}, m={p['m']}, n={p['n']}, Rw={p['rw']})", fontsize=9); ax2.grid(alpha=.3, which='both')
        fig.tight_layout(); self.canv['Crossplots'].draw()

    def summary(self):
        p, d, r = self.used, self.w['depth'], self.res
        L = ['INTERPRETATION SUMMARY', f"Well: {self.w['well']}", '', 'PARAMETERS USED (stored with results):']
        L += [f'  {k}: {v}' for k, v in p.items()]
        L += ['', 'ZONE STATISTICS (calculated; mean [min-max]):']
        for line in self.zones.get('1.0', 'end').strip().splitlines():
            try: n, t, b = [x.strip() for x in line.split(',')]; z = P.zone_stats(d, r, float(t), float(b))
            except Exception: L.append(f'  (could not parse zone line: {line})'); continue
            f = lambda s: f'{s[0]:.3f} [{s[1]:.3f}-{s[2]:.3f}]'
            L += [f"  {n}: {t}-{b} {self.w['depth_unit']}", f"    Gross {z['gross']:.1f}  Net pay {z['net']:.1f}  NTG {z['ntg']:.2f}",
                  f"    Vsh {f(z['vsh'])}", f"    PHI {f(z['phi'])}", f"    Sw  {f(z['sw'])}"]
        L += ['', 'NOTE: Permeability not estimated. Pay flags are criteria-based screening, not a definitive result.', DISCLAIMER]
        self.settext('Zones / Summary', '\n'.join(L))

    # ---------- export ----------
    def need(self):
        if self.res is None: messagebox.showinfo('Nothing to export', 'Run CALCULATE first.'); return False
        return True

    def exp_csv(self):
        if not self.need(): return
        p = filedialog.asksaveasfilename(defaultextension='.csv')
        if not p: return
        names = self.w['names']; r = self.res
        with open(p, 'w') as f:
            f.write('# ' + '; '.join(f'{k}={v}' for k, v in self.used.items()) + '\n')
            f.write('DEPTH,' + ','.join(names) + ',VSH,PHI,SW,SH,PAY\n')
            for i, dep in enumerate(self.w['depth']):
                f.write(f'{dep},' + ','.join(f"{self.w['data'][n][i]}" for n in names) +
                        f",{r['VSH'][i]},{r['PHI'][i]},{r['SW'][i]},{r['SH'][i]},{int(r['PAY'][i])}\n")

    def exp_png(self):
        tab = self.nb.tab(self.nb.select(), 'text')
        if tab not in self.figs: messagebox.showinfo('PNG', 'Select the Log Viewer or Crossplots tab.'); return
        p = filedialog.asksaveasfilename(defaultextension='.png')
        if p: self.figs[tab].savefig(p, dpi=200)

    def exp_pdf(self):
        if not self.need(): return
        p = filedialog.asksaveasfilename(defaultextension='.pdf')
        if not p: return
        with PdfPages(p) as pdf:
            tf = Figure(figsize=(8.5, 11)); tf.text(.05, .97, self.txt['QC'].get('1.0', 'end') + '\n' +
                                                     self.txt['Zones / Summary'].get('1.0', 'end'), va='top', fontsize=6, family='monospace')
            pdf.savefig(tf)
            for f in self.figs.values(): pdf.savefig(f)


if __name__ == '__main__':
    App().mainloop()

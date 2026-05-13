# Optional Gnuplot Panel Scripts

These scripts plot individual panels using the final data in `../outputs/`.
The main PI-ready gnuplot figure is still:

```bash
cd ..
gnuplot fig11_final_hex.gnu
```

Use these only if you want to inspect or export one panel at a time:

```bash
cd ..
gnuplot gnuplot_panels/fig11_panel_a.gnu
gnuplot gnuplot_panels/fig11_panel_b.gnu
gnuplot gnuplot_panels/fig11_panel_c.gnu
```

The scripts assume they are launched from the main `triangular/` folder so the
relative paths to `outputs/` resolve correctly.

import { Component, computed, inject } from '@angular/core';
import { RouterOutlet } from '@angular/router';
import { ButtonModule } from 'primeng/button';
import { CardModule } from 'primeng/card';
import { ToolbarModule } from 'primeng/toolbar';

import { ThemeService } from './core/theme.service';

@Component({
  selector: 'app-root',
  imports: [RouterOutlet, ToolbarModule, ButtonModule, CardModule],
  templateUrl: './app.html',
  styleUrl: './app.scss',
})
export class App {
  private readonly theme = inject(ThemeService);

  protected readonly darkMode = this.theme.darkMode;

  protected readonly themeIcon = computed(() => (this.darkMode() ? 'pi pi-sun' : 'pi pi-moon'));

  protected toggleTheme(): void {
    this.theme.toggle();
  }
}

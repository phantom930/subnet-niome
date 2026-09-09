import { Injectable, effect, signal } from '@angular/core';

const STORAGE_KEY = 'niome-dashboard-theme';
const DARK_CLASS = 'app-dark';

/**
 * Tracks light/dark mode and toggles the class that PrimeNG's
 * `darkModeSelector` watches. The choice is remembered per browser.
 */
@Injectable({ providedIn: 'root' })
export class ThemeService {
  readonly darkMode = signal<boolean>(this.readInitialPreference());

  constructor() {
    effect(() => {
      const dark = this.darkMode();
      document.documentElement.classList.toggle(DARK_CLASS, dark);
      try {
        localStorage.setItem(STORAGE_KEY, dark ? 'dark' : 'light');
      } catch {
        // Private browsing or blocked storage: keep the in-memory value only.
      }
    });
  }

  toggle(): void {
    this.darkMode.update((dark) => !dark);
  }

  private readInitialPreference(): boolean {
    try {
      const stored = localStorage.getItem(STORAGE_KEY);
      if (stored === 'dark') return true;
      if (stored === 'light') return false;
    } catch {
      // Fall through to the OS preference.
    }
    return window.matchMedia?.('(prefers-color-scheme: dark)').matches ?? false;
  }
}
